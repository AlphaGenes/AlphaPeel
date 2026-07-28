"""Main peeling functions."""

import concurrent.futures

from numba import jit, float32
import numpy as np


_LOCUS_THREAD_COUNT = 1
_DEFER_SEGREGATION_COLLAPSE = False


def set_locus_thread_count(n_threads):
    """Set the number of Python threads used for manual locus splitting."""

    global _LOCUS_THREAD_COUNT
    _LOCUS_THREAD_COUNT = n_threads


def set_defer_segregation_collapse(enabled):
    """Set whether peel-down estimates segregation without collapsing it."""

    global _DEFER_SEGREGATION_COLLAPSE
    _DEFER_SEGREGATION_COLLAPSE = enabled


def peel_down(family, peeling_info, single_locus_mode):
    """Peel information down from parents to offspring."""

    n_threads = _LOCUS_THREAD_COUNT
    if n_threads <= 1:
        _peel_down_serial(
            family, peeling_info, single_locus_mode, _DEFER_SEGREGATION_COLLAPSE
        )
        return

    _peel_down_threaded_parent_setup(family, peeling_info, single_locus_mode, n_threads)


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e16": float32, "e1e": float32},
)
def _peel_down_serial(family, peeling_info, single_locus_mode, defer_collapse):
    """Peel information down from parents to offspring.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param single_locus_mode: A flag to indicate the mode of peeling.
        `False` if using multi-locus peeling, and `True` if using single-locus peeling.
    :type single_locus_mode: bool
    :return: None. The function modifies the peeling_info object in place.
    """
    e = 0.000001
    e1e = 1 - e
    e4 = e / 4
    e16 = e / 16

    n_loci = peeling_info.n_loci
    n_offspring = len(family.offspring)
    workspace = create_peel_down_workspace(n_offspring, n_loci)
    # workspace[2]: prob_sire, workspace[3]: prob_dam.
    prob_sire = workspace[2]
    prob_dam = workspace[3]

    # Build parent genotype probabilities excluding this family's current contribution.
    setup_parent_genotype_probs(
        family,
        peeling_info,
        (prob_sire, prob_dam),
        n_loci,
    )

    # Equivalent einsum: joint_parents = np.einsum("ai, bi -> abi", prob_sire, prob_dam)
    joint_parents = get_joint_parents(prob_sire, prob_dam, n_loci)
    smooth_16_by_locus(joint_parents, e1e, e16, n_loci)

    accumulate_children_for_peel_down(
        family, peeling_info, workspace, (e1e, e4), single_locus_mode
    )

    # Combine all child projections on the log scale, then subtract each child
    # to produce child-specific parent estimates.
    add_joint_parents_and_all_to_minus(
        workspace[7], joint_parents, workspace[1], n_offspring, n_loci
    )
    normalize_parent_estimates_for_peel_down(workspace, n_offspring, n_loci)
    update_child_anteriors_for_peel_down(family, peeling_info, workspace, n_loci)

    if not single_locus_mode:
        update_child_segregation_for_peel_down(
            family, peeling_info, workspace, (e1e, e4), defer_collapse
        )


def _peel_down_threaded_parent_setup(
    family, peeling_info, single_locus_mode, n_threads
):
    """Peel down with parent/joint setup split across locus chunks."""

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4
    e16 = e / 16

    n_loci = peeling_info.n_loci
    n_offspring = len(family.offspring)
    workspace = create_peel_down_workspace(n_offspring, n_loci)
    prob_sire = workspace[2]
    prob_dam = workspace[3]
    joint_parents = np.empty((4, 4, n_loci), dtype=np.float32)

    _setup_parent_probs_and_joint_by_locus_threads(
        family,
        peeling_info,
        prob_sire,
        prob_dam,
        joint_parents,
        e1e,
        e16,
        n_loci,
        n_threads,
    )

    _peel_down_after_parent_setup_by_locus_threads(
        family,
        peeling_info,
        workspace,
        joint_parents,
        (e1e, e4),
        single_locus_mode,
        _DEFER_SEGREGATION_COLLAPSE,
        n_loci,
        n_offspring,
        n_threads,
    )

    if not single_locus_mode and not _DEFER_SEGREGATION_COLLAPSE:
        collapse_estimated_child_segregations_for_peel_down(
            family, peeling_info, workspace, (e1e, e4)
        )


def _setup_parent_probs_and_joint_by_locus_threads(
    family,
    peeling_info,
    prob_sire,
    prob_dam,
    joint_parents,
    joint_scale,
    joint_floor,
    n_loci,
    n_threads,
):
    """Split the parent/joint setup by locus and run chunks on Python threads."""

    n_workers = min(n_threads, n_loci)
    chunk_size = (n_loci + n_workers - 1) // n_workers
    futures = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for start in range(0, n_loci, chunk_size):
            stop = min(start + chunk_size, n_loci)
            futures.append(
                executor.submit(
                    setup_parent_probs_and_joint_slice,
                    peeling_info.posterior[family.sire, :, start:stop],
                    peeling_info.posterior[family.dam, :, start:stop],
                    peeling_info.posterior_sire_contribution[family.idn, :, start:stop],
                    peeling_info.posterior_dam_contribution[family.idn, :, start:stop],
                    peeling_info.anterior[family.sire, :, start:stop],
                    peeling_info.anterior[family.dam, :, start:stop],
                    peeling_info.penetrance[family.sire, :, start:stop],
                    peeling_info.penetrance[family.dam, :, start:stop],
                    prob_sire[:, start:stop],
                    prob_dam[:, start:stop],
                    joint_parents[:, :, start:stop],
                    joint_scale,
                    joint_floor,
                    stop - start,
                )
            )

        for future in futures:
            future.result()


def _peel_down_after_parent_setup_by_locus_threads(
    family,
    peeling_info,
    workspace,
    joint_parents,
    smoothing,
    single_locus_mode,
    defer_collapse,
    n_loci,
    n_offspring,
    n_threads,
):
    """Split the post-parent peel-down calculations by locus."""

    n_workers = min(n_threads, n_loci)
    chunk_size = (n_loci + n_workers - 1) // n_workers
    futures = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for start in range(0, n_loci, chunk_size):
            stop = min(start + chunk_size, n_loci)
            futures.append(
                executor.submit(
                    peel_down_after_parent_setup_slice,
                    family,
                    peeling_info,
                    workspace,
                    joint_parents[:, :, start:stop],
                    smoothing[0],
                    smoothing[1],
                    single_locus_mode,
                    start,
                    stop,
                    n_offspring,
                )
            )

        for future in futures:
            future.result()


@jit(nopython=True, nogil=True)
def create_peel_down_workspace(n_offspring, n_loci):
    """Create reusable arrays for the peel-down pass."""

    return (
        # workspace[0]: current child's projection onto joint parent genotypes.
        np.full((4, 4, n_loci), 0, dtype=np.float32),
        # workspace[1]: all children's projections accumulated on joint parents.
        np.full((4, 4, n_loci), 0, dtype=np.float32),
        # workspace[2]: sire genotype probabilities with this family removed.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[3]: dam genotype probabilities with this family removed.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[4]: current child's posterior and penetrance by genotype state.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[5]: current child's segregation probabilities.
        np.full((4, n_loci), 1, dtype=np.float32),
        # workspace[6]: per-child inheritance tensors for anterior updates.
        np.full((n_offspring, 4, 4, 4, n_loci), 0, dtype=np.float32),
        # workspace[7]: joint parent estimates with one child contribution removed.
        np.full((n_offspring, 4, 4, n_loci), 1, dtype=np.float32),
        # workspace[8]: stored child genotype values for segregation updates.
        np.full((n_offspring, 4, n_loci), 0, dtype=np.float32),
        # workspace[9]: forward-pass buffer for collapsed segregation updates.
        np.full((4, n_loci), 1, dtype=np.float32),
    )


@jit(nopython=True, nogil=True)
def accumulate_children_for_peel_down(
    family, peeling_info, workspace, smoothing, single_locus_mode
):
    """Project children onto joint parent genotypes for peel-down."""

    n_loci = peeling_info.n_loci
    # workspace[0]: child_to_parents_current, [1]: all_to_parents,
    # [4]: child_values_current, [5]: current_seg, [6]: child_seg_tensor,
    # [7]: parents_minus_child, [8]: child_values_tensor.
    for index, child in enumerate(family.offspring):
        child_segs = workspace[6][index, :, :, :, :]
        if peeling_info.is_x_chr:
            project_child_for_peel_x_chr(
                child,
                peeling_info,
                (workspace[4], workspace[5], child_segs, workspace[0]),
                n_loci,
                smoothing,
            )
        else:
            project_child_for_peel_autosome(
                child,
                peeling_info,
                (workspace[4], workspace[5], child_segs, workspace[0]),
                n_loci,
                smoothing,
            )
        if not single_locus_mode:
            workspace[8][index, :, :] = workspace[4]

        add_log_child_to_parents_and_subtract_current(
            workspace[0], workspace[1], workspace[7][index, :, :, :], n_loci
        )


@jit(nopython=True, nogil=True)
def normalize_parent_estimates_for_peel_down(workspace, n_offspring, n_loci):
    """Convert peel-down parent estimates from log scale to probabilities."""

    # workspace[1]: all_to_parents, workspace[7]: parents_minus_child.
    workspace[1][:, :, :] = exp_norm_2d(workspace[1], n_loci)
    for i in range(n_offspring):
        workspace[7][i, :, :, :] = exp_norm_2d(workspace[7][i, :, :, :], n_loci)


@jit(nopython=True, nogil=True)
def update_child_anteriors_for_peel_down(family, peeling_info, workspace, n_loci):
    """Update child anterior genotype probabilities after peel-down."""

    # workspace[6]: child_seg_tensor, workspace[7]: parents_minus_child.
    for index, child in enumerate(family.offspring):
        child_anterior = peeling_info.anterior[child, :, :]
        # Equivalent einsum: anterior[child, :, :] =
        # np.einsum("abci, abi -> ci", child_seg_tensor[i, :, :, :, :],
        # parents_minus_child[i, :, :, :])
        project_parent_genotypes(
            workspace[6][index, :, :, :, :],
            workspace[7][index, :, :, :],
            child_anterior,
            n_loci,
        )
        normalize_4_by_locus(child_anterior, n_loci)


@jit(nopython=True, nogil=True)
def normalize_parent_estimates_slice(workspace, start, stop, n_offspring):
    """Convert sliced peel-down parent estimates from log scale to probabilities."""

    n_loci = stop - start
    workspace[1][:, :, start:stop] = exp_norm_2d(workspace[1][:, :, start:stop], n_loci)
    for i in range(n_offspring):
        workspace[7][i, :, :, start:stop] = exp_norm_2d(
            workspace[7][i, :, :, start:stop], n_loci
        )


@jit(nopython=True, nogil=True)
def update_child_anteriors_for_peel_down_slice(
    family, peeling_info, workspace, start, stop
):
    """Update child anterior genotype probabilities for one locus slice."""

    n_loci = stop - start
    for index, child in enumerate(family.offspring):
        project_parent_genotypes(
            workspace[6][index, :, :, :, start:stop],
            workspace[7][index, :, :, start:stop],
            peeling_info.anterior[child, :, start:stop],
            n_loci,
        )
        normalize_4_by_locus(peeling_info.anterior[child, :, start:stop], n_loci)


@jit(nopython=True, nogil=True)
def update_child_segregation_for_peel_down(
    family, peeling_info, workspace, smoothing, defer_collapse
):
    """Update child segregation probabilities after peel-down."""

    for index, child in enumerate(family.offspring):
        update_one_child_segregation_for_peel_down(
            child, peeling_info, workspace, index, smoothing, defer_collapse
        )


@jit(nopython=True, nogil=True)
def collapse_estimated_child_segregations_for_peel_down(
    family, peeling_info, workspace, smoothing
):
    """Collapse already-estimated child segregations after threaded locus work."""

    for child in family.offspring:
        collapse_one_child_segregation_for_peel_down(child, peeling_info, workspace)
        smooth_4_by_locus(
            peeling_info.segregation[child, :, :],
            smoothing[0],
            smoothing[1],
            peeling_info.n_loci,
        )


@jit(nopython=True, nogil=True)
def update_one_child_segregation_for_peel_down(
    child, peeling_info, workspace, child_index, smoothing, defer_collapse
):
    """Update one child's segregation probabilities after peel-down."""

    estimate_one_child_segregation_for_peel_down(
        child, peeling_info, workspace, child_index
    )
    if defer_collapse:
        return

    collapse_one_child_segregation_for_peel_down(child, peeling_info, workspace)
    smooth_4_by_locus(
        peeling_info.segregation[child, :, :],
        smoothing[0],
        smoothing[1],
        peeling_info.n_loci,
    )


@jit(nopython=True, nogil=True)
def estimate_one_child_segregation_for_peel_down(
    child, peeling_info, workspace, child_index
):
    """Estimate one child's uncollapsed segregation probabilities."""

    if peeling_info.is_x_chr and peeling_info.sex[child] == 0:  # 0=male, 1=female.
        segregation_tensor = peeling_info.segregation_tensor_xy
        segregation_tensor_norm = peeling_info.segregation_tensor_xy_norm
    elif peeling_info.is_x_chr:
        segregation_tensor = peeling_info.segregation_tensor_xx
        segregation_tensor_norm = peeling_info.segregation_tensor_xx_norm
    else:
        segregation_tensor = peeling_info.segregation_tensor
        segregation_tensor_norm = peeling_info.segregation_tensor_norm

    # Equivalent einsum: segregation[child, :, :] =
    # np.einsum("abcd, abi, ci -> di", segregation_tensor,
    # parents_minus_child, child_values)
    # workspace[7]: parents_minus_child, workspace[8]: child_values_tensor.
    estimate_segregation_with_norm(
        (segregation_tensor, segregation_tensor_norm),
        (workspace[7][child_index, :, :, :], workspace[8][child_index, :, :]),
        peeling_info.segregation[child, :, :],
        peeling_info.n_loci,
    )


@jit(nopython=True, nogil=True)
def estimate_child_segregations_for_peel_down_slice(
    family, peeling_info, workspace, start, stop
):
    """Estimate child segregation probabilities for one locus slice."""

    n_loci = stop - start
    for index, child in enumerate(family.offspring):
        if peeling_info.is_x_chr and peeling_info.sex[child] == 0:  # 0=male, 1=female.
            segregation_tensor = peeling_info.segregation_tensor_xy
            segregation_tensor_norm = peeling_info.segregation_tensor_xy_norm
        elif peeling_info.is_x_chr:
            segregation_tensor = peeling_info.segregation_tensor_xx
            segregation_tensor_norm = peeling_info.segregation_tensor_xx_norm
        else:
            segregation_tensor = peeling_info.segregation_tensor
            segregation_tensor_norm = peeling_info.segregation_tensor_norm

        estimate_segregation_with_norm(
            (segregation_tensor, segregation_tensor_norm),
            (
                workspace[7][index, :, :, start:stop],
                workspace[8][index, :, start:stop],
            ),
            peeling_info.segregation[child, :, start:stop],
            n_loci,
        )


@jit(nopython=True, nogil=True)
def collapse_one_child_segregation_for_peel_down(child, peeling_info, workspace):
    """Collapse one child's full-chromosome segregation probabilities."""

    # workspace[9]: forward_seg
    collapse_segregation_in_place(
        peeling_info.segregation[child, :, :],
        peeling_info.transmission_rate,
        workspace[9],
        peeling_info.n_loci,
    )


@jit(nopython=True, nogil=True)
def collapse_child_segregations_in_place(
    segregation, transmission, child_ids, scale, floor, n_loci
):
    """Collapse full-chromosome segregation probabilities for selected children."""

    forward = np.full((4, n_loci), 1, dtype=np.float32)
    for child in child_ids:
        collapse_segregation_in_place(
            segregation[child, :, :],
            transmission,
            forward,
            n_loci,
        )
        smooth_4_by_locus(segregation[child, :, :], scale, floor, n_loci)


@jit(nopython=True, nogil=True)
def peel_down_after_parent_setup_slice(
    family,
    peeling_info,
    workspace,
    joint_parents,
    scale,
    floor,
    single_locus_mode,
    start,
    stop,
    n_offspring,
):
    """Run post-parent peel-down calculations for one locus slice."""

    n_loci = stop - start
    child_to_parents = workspace[0][:, :, start:stop]
    all_to_parents = workspace[1][:, :, start:stop]
    child_values = workspace[4][:, start:stop]
    current_seg = workspace[5][:, start:stop]

    for index, child in enumerate(family.offspring):
        child_segs = workspace[6][index, :, :, :, start:stop]
        if peeling_info.is_x_chr and peeling_info.sex[child] == 0:
            segregation_tensor = peeling_info.segregation_tensor_xy
        elif peeling_info.is_x_chr:
            segregation_tensor = peeling_info.segregation_tensor_xx
        else:
            segregation_tensor = peeling_info.segregation_tensor

        project_child_for_peel_slice(
            peeling_info.posterior[child, :, start:stop],
            peeling_info.penetrance[child, :, start:stop],
            peeling_info.segregation[child, :, start:stop],
            segregation_tensor,
            child_values,
            current_seg,
            child_segs,
            child_to_parents,
            scale,
            floor,
            n_loci,
        )
        if not single_locus_mode:
            workspace[8][index, :, start:stop] = child_values

        add_log_child_to_parents_and_subtract_current(
            child_to_parents,
            all_to_parents,
            workspace[7][index, :, :, start:stop],
            n_loci,
        )

    add_joint_parents_and_all_to_minus(
        workspace[7][:, :, :, start:stop],
        joint_parents,
        all_to_parents,
        n_offspring,
        n_loci,
    )
    normalize_parent_estimates_slice(workspace, start, stop, n_offspring)
    update_child_anteriors_for_peel_down_slice(
        family, peeling_info, workspace, start, stop
    )

    if not single_locus_mode:
        estimate_child_segregations_for_peel_down_slice(
            family, peeling_info, workspace, start, stop
        )


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e1e": float32},
)
def _peel_up_serial(family, peeling_info):
    """Peel information up from offspring to parents.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function modifies the peeling_info object in place.
    """
    e = 0.000001
    e1e = 1 - e
    e4 = e / 4

    n_loci = peeling_info.n_loci
    workspace = create_peel_up_workspace(n_loci)

    # workspace[2]: prob_sire, workspace[3]: prob_dam.
    setup_parent_genotype_probs(
        family,
        peeling_info,
        (workspace[2], workspace[3]),
        n_loci,
    )
    smooth_4_by_locus(workspace[2], e1e, e4, n_loci)
    smooth_4_by_locus(workspace[3], e1e, e4, n_loci)
    accumulate_children_for_peel_up(family, peeling_info, workspace, (e1e, e4))
    update_parent_posteriors_for_peel_up(family, peeling_info, workspace, (e1e, e4))


def peel_up(family, peeling_info):
    """Peel information up from offspring to parents."""

    n_threads = _LOCUS_THREAD_COUNT
    if n_threads <= 1:
        _peel_up_serial(family, peeling_info)
        return

    _peel_up_threaded_by_locus(family, peeling_info, n_threads)


def _peel_up_threaded_by_locus(family, peeling_info, n_threads):
    """Peel up with each Python thread owning a locus chunk."""

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4

    n_loci = peeling_info.n_loci
    n_workers = min(n_threads, n_loci)
    chunk_size = (n_loci + n_workers - 1) // n_workers
    workspace = create_peel_up_workspace(n_loci)
    futures = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for start in range(0, n_loci, chunk_size):
            stop = min(start + chunk_size, n_loci)
            futures.append(
                executor.submit(
                    _peel_up_locus_slice,
                    family,
                    peeling_info,
                    workspace,
                    e1e,
                    e4,
                    start,
                    stop,
                )
            )

        for future in futures:
            future.result()


def _peel_up_locus_slice(
    family,
    peeling_info,
    workspace,
    scale,
    floor,
    start,
    stop,
):
    """Run the complete peel-up calculation for one locus slice."""

    n_loci = stop - start
    prob_sire = workspace[2][:, start:stop]
    prob_dam = workspace[3][:, start:stop]
    child_values = workspace[4][:, start:stop]
    current_seg = workspace[5][:, start:stop]
    child_segs = workspace[6][0, :, :, :, start:stop]
    child_to_parents = workspace[0][:, :, start:stop]
    all_to_parents = workspace[1][:, :, start:stop]

    setup_parent_probs_slice(
        peeling_info.posterior[family.sire, :, start:stop],
        peeling_info.posterior[family.dam, :, start:stop],
        peeling_info.posterior_sire_contribution[family.idn, :, start:stop],
        peeling_info.posterior_dam_contribution[family.idn, :, start:stop],
        peeling_info.anterior[family.sire, :, start:stop],
        peeling_info.anterior[family.dam, :, start:stop],
        peeling_info.penetrance[family.sire, :, start:stop],
        peeling_info.penetrance[family.dam, :, start:stop],
        prob_sire,
        prob_dam,
        n_loci,
    )
    smooth_4_by_locus(prob_sire, scale, floor, n_loci)
    smooth_4_by_locus(prob_dam, scale, floor, n_loci)

    for child in family.offspring:
        if peeling_info.is_x_chr and peeling_info.sex[child] == 0:
            segregation_tensor = peeling_info.segregation_tensor_xy
        elif peeling_info.is_x_chr:
            segregation_tensor = peeling_info.segregation_tensor_xx
        else:
            segregation_tensor = peeling_info.segregation_tensor

        project_child_for_peel_slice(
            peeling_info.posterior[child, :, start:stop],
            peeling_info.penetrance[child, :, start:stop],
            peeling_info.segregation[child, :, start:stop],
            segregation_tensor,
            child_values,
            current_seg,
            child_segs,
            child_to_parents,
            scale,
            floor,
            n_loci,
        )
        add_log_child_to_parents(child_to_parents, all_to_parents, n_loci)

    update_parent_posteriors_for_peel_up_slice(
        all_to_parents,
        prob_sire,
        prob_dam,
        peeling_info.posterior_sire_contribution[family.idn, :, start:stop],
        peeling_info.posterior_dam_contribution[family.idn, :, start:stop],
        scale,
        floor,
        n_loci,
    )


@jit(nopython=True, nogil=True)
def create_peel_up_workspace(n_loci):
    """Create reusable arrays for the peel-up pass."""

    return (
        # workspace[0]: current child's projection onto joint parent genotypes.
        np.full((4, 4, n_loci), 0, dtype=np.float32),
        # workspace[1]: all children's projections accumulated on joint parents.
        np.full((4, 4, n_loci), 0, dtype=np.float32),
        # workspace[2]: sire genotype probabilities with this family removed.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[3]: dam genotype probabilities with this family removed.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[4]: current child's posterior and penetrance by genotype state.
        np.full((4, n_loci), 0, dtype=np.float32),
        # workspace[5]: current child's segregation probabilities.
        np.full((4, n_loci), 1, dtype=np.float32),
        # workspace[6]: single-child inheritance tensor reused for each offspring.
        np.full((1, 4, 4, 4, n_loci), 0, dtype=np.float32),
    )


@jit(nopython=True, nogil=True)
def accumulate_children_for_peel_up(family, peeling_info, workspace, smoothing):
    """Project children onto joint parent genotypes for peel-up."""

    n_loci = peeling_info.n_loci
    # workspace[0]: child_to_parents_current, [1]: all_to_parents,
    # [4]: child_values_current, [5]: current_seg,
    # [6]: forward_seg.
    child_workspace = (
        workspace[4],
        workspace[5],
        workspace[6][0, :, :, :, :],
        workspace[0],
    )
    for child in family.offspring:
        if peeling_info.is_x_chr:
            project_child_for_peel_x_chr(
                child, peeling_info, child_workspace, n_loci, smoothing
            )
        else:
            project_child_for_peel_autosome(
                child, peeling_info, child_workspace, n_loci, smoothing
            )
        add_log_child_to_parents(workspace[0], workspace[1], n_loci)


@jit(nopython=True, nogil=True)
def update_parent_posteriors_for_peel_up(family, peeling_info, workspace, smoothing):
    """Update parent posterior contributions after peel-up."""

    n_loci = peeling_info.n_loci
    # workspace[1]: all_to_parents.
    all_to_parents = exp_norm_2d(workspace[1], n_loci)

    sire_posterior = peeling_info.posterior_sire_contribution[family.idn, :, :]
    # workspace[3]: prob_dam.
    combine_and_reduce_axis1(all_to_parents, workspace[3], sire_posterior, n_loci)
    normalize_4_by_locus(sire_posterior, n_loci)
    smooth_4_by_locus(sire_posterior, smoothing[0], smoothing[1], n_loci)

    dam_posterior = peeling_info.posterior_dam_contribution[family.idn, :, :]
    # workspace[2]: prob_sire.
    combine_and_reduce_axis0(all_to_parents, workspace[2], dam_posterior, n_loci)
    normalize_4_by_locus(dam_posterior, n_loci)
    smooth_4_by_locus(dam_posterior, smoothing[0], smoothing[1], n_loci)


@jit(nopython=True, nogil=True)
def setup_parent_probs_slice(
    posterior_sire,
    posterior_dam,
    posterior_sire_contribution,
    posterior_dam_contribution,
    anterior_sire,
    anterior_dam,
    penetrance_sire,
    penetrance_dam,
    prob_sire,
    prob_dam,
    n_loci,
):
    """Build normalized parent probabilities for a locus slice."""

    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] = np.log(posterior_sire[state, locus]) - np.log(
                posterior_sire_contribution[state, locus]
            )
            prob_dam[state, locus] = np.log(posterior_dam[state, locus]) - np.log(
                posterior_dam_contribution[state, locus]
            )

    exp_norm_1d_in_place(prob_sire, n_loci)
    exp_norm_1d_in_place(prob_dam, n_loci)

    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] *= (
                anterior_sire[state, locus] * penetrance_sire[state, locus]
            )
            prob_dam[state, locus] *= (
                anterior_dam[state, locus] * penetrance_dam[state, locus]
            )

    normalize_4_by_locus(prob_sire, n_loci)
    normalize_4_by_locus(prob_dam, n_loci)


@jit(nopython=True, nogil=True)
def setup_parent_probs_and_joint_slice(
    posterior_sire,
    posterior_dam,
    posterior_sire_contribution,
    posterior_dam_contribution,
    anterior_sire,
    anterior_dam,
    penetrance_sire,
    penetrance_dam,
    prob_sire,
    prob_dam,
    joint_parents,
    joint_scale,
    joint_floor,
    n_loci,
):
    """Build parent probabilities and smoothed joint estimates for a locus slice."""

    setup_parent_probs_slice(
        posterior_sire,
        posterior_dam,
        posterior_sire_contribution,
        posterior_dam_contribution,
        anterior_sire,
        anterior_dam,
        penetrance_sire,
        penetrance_dam,
        prob_sire,
        prob_dam,
        n_loci,
    )

    for geno_sire in range(4):
        for geno_dam in range(4):
            for locus in range(n_loci):
                joint_parents[geno_sire, geno_dam, locus] = (
                    prob_sire[geno_sire, locus]
                    * prob_dam[geno_dam, locus]
                    * joint_scale
                    + joint_floor
                )


@jit(nopython=True, nogil=True)
def project_child_for_peel_slice(
    posterior_child,
    penetrance_child,
    segregation_child,
    segregation_tensor,
    child_values,
    current_seg,
    child_segs,
    child_to_parents,
    scale,
    floor,
    n_loci,
):
    """Project one child's sliced loci onto parental genotypes."""

    for state in range(4):
        for locus in range(n_loci):
            child_values[state, locus] = (
                posterior_child[state, locus] * penetrance_child[state, locus]
            )
            current_seg[state, locus] = segregation_child[state, locus]

    normalize_4_by_locus(child_values, n_loci)
    smooth_4_by_locus(child_values, scale, floor, n_loci)
    normalize_4_by_locus(current_seg, n_loci)

    create_child_segs(segregation_tensor, current_seg, child_segs, n_loci)
    project_child_genotypes(child_segs, child_values, child_to_parents, n_loci)


@jit(nopython=True, nogil=True)
def update_parent_posteriors_for_peel_up_slice(
    all_to_parents,
    prob_sire,
    prob_dam,
    sire_posterior,
    dam_posterior,
    scale,
    floor,
    n_loci,
):
    """Update parent posterior contributions for a locus slice."""

    all_to_parents = exp_norm_2d(all_to_parents, n_loci)

    combine_and_reduce_axis1(all_to_parents, prob_dam, sire_posterior, n_loci)
    normalize_4_by_locus(sire_posterior, n_loci)
    smooth_4_by_locus(sire_posterior, scale, floor, n_loci)

    combine_and_reduce_axis0(all_to_parents, prob_sire, dam_posterior, n_loci)
    normalize_4_by_locus(dam_posterior, n_loci)
    smooth_4_by_locus(dam_posterior, scale, floor, n_loci)


@jit(nopython=True, nogil=True)
def setup_parent_genotype_probs(
    family,
    peeling_info,
    parent_probs,
    n_loci,
):
    """Build the normalized sire and dam genotype probabilities for a family."""

    prob_sire, prob_dam = parent_probs
    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] = np.log(
                peeling_info.posterior[family.sire, state, locus]
            ) - np.log(
                peeling_info.posterior_sire_contribution[family.idn, state, locus]
            )
            prob_dam[state, locus] = np.log(
                peeling_info.posterior[family.dam, state, locus]
            ) - np.log(
                peeling_info.posterior_dam_contribution[family.idn, state, locus]
            )

    exp_norm_1d_in_place(prob_sire, n_loci)
    exp_norm_1d_in_place(prob_dam, n_loci)

    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] *= (
                peeling_info.anterior[family.sire, state, locus]
                * peeling_info.penetrance[family.sire, state, locus]
            )
            prob_dam[state, locus] *= (
                peeling_info.anterior[family.dam, state, locus]
                * peeling_info.penetrance[family.dam, state, locus]
            )

    normalize_4_by_locus(prob_sire, n_loci)
    normalize_4_by_locus(prob_dam, n_loci)


@jit(
    nopython=True,
    nogil=True,
    locals={"e4": float32, "e1e": float32},
)
def project_child_for_peel_autosome(
    child,
    peeling_info,
    child_workspace,
    n_loci,
    smoothing,
):
    """Project one autosomal child onto parental genotypes."""

    child_values, current_seg, child_segs, child_to_parents = child_workspace
    for state in range(4):
        for locus in range(n_loci):
            child_values[state, locus] = (
                peeling_info.posterior[child, state, locus]
                * peeling_info.penetrance[child, state, locus]
            )
            current_seg[state, locus] = peeling_info.segregation[child, state, locus]

    normalize_4_by_locus(child_values, n_loci)
    smooth_4_by_locus(child_values, smoothing[0], smoothing[1], n_loci)
    normalize_4_by_locus(current_seg, n_loci)

    create_child_segs(peeling_info.segregation_tensor, current_seg, child_segs, n_loci)
    project_child_genotypes(child_segs, child_values, child_to_parents, n_loci)


@jit(
    nopython=True,
    nogil=True,
    locals={"e4": float32, "e1e": float32},
)
def project_child_for_peel_x_chr(
    child,
    peeling_info,
    child_workspace,
    n_loci,
    smoothing,
):
    """Project one X-chromosome child onto parental genotypes."""

    child_values, current_seg, child_segs, child_to_parents = child_workspace
    for state in range(4):
        for locus in range(n_loci):
            child_values[state, locus] = (
                peeling_info.posterior[child, state, locus]
                * peeling_info.penetrance[child, state, locus]
            )
            current_seg[state, locus] = peeling_info.segregation[child, state, locus]

    normalize_4_by_locus(child_values, n_loci)
    smooth_4_by_locus(child_values, smoothing[0], smoothing[1], n_loci)
    normalize_4_by_locus(current_seg, n_loci)

    if peeling_info.sex[child] == 0:  # 0=male, 1=female.
        segregation_tensor = peeling_info.segregation_tensor_xy
    else:
        segregation_tensor = peeling_info.segregation_tensor_xx

    create_child_segs(segregation_tensor, current_seg, child_segs, n_loci)
    project_child_genotypes(child_segs, child_values, child_to_parents, n_loci)


@jit(nopython=True, nogil=True)
def get_joint_parents(prob_sire, prob_dam, n_loci):
    """Creates the joint parental genotypes based on the probabilities for each parent.

    :param prob_sire: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type prob_sire: 2D numpy array of float32 with size 4 x n_loci
    :param prob_dam: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type prob_dam: 2D numpy array of float32 with size 4 x n_loci
    :return: the probabilities of all the combinations of the sire's genotype and the dam's genotype
        with information from previous peeling cycle (P(p, m))
    :rtype: 3D numpy array of float32 with size 4 x 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("ai, bi -> abi", prob_sire, prob_dam)
    output = np.full(shape=(4, 4, n_loci), fill_value=0, dtype=np.float32)
    for geno_sire in range(4):
        for geno_dam in range(4):
            for locus in range(n_loci):
                output[geno_sire, geno_dam, locus] = (
                    prob_sire[geno_sire, locus] * prob_dam[geno_dam, locus]
                )
    return output


@jit(nopython=True, nogil=True)
def create_child_segs(segregation_tensor, current_seg, output, n_loci):
    """Creates the child-specific segregation tensor using the child's current segregation estimate.

    :param segregation_tensor: the probability of each combination of
        the sire's genotype, the dam's genotype, child's genotype and
        segregation without any other information (P(p, m, allele, seg))
    :type segregation_tensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param current_seg: The probability of each segregation of each locus of the child
        with information from previous peeling cycle (P(seg))
    :type current_seg: 2D numpy array of float32 with size 4 x n_loci
    :param output: the probability of each combination of the sire's genotype,
        the dam's genotype and the child's genotype of each locus with information
        from previous peeling cycle (P(p, m, allele))
    :type output: 4D numpy array of float32 with size 4 x 4 x 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abcd, di -> abci", segregation_tensor, current_seg)
    for geno_sire in range(4):
        for geno_dam in range(4):
            for geno_child in range(4):
                for locus in range(n_loci):
                    total = 0
                    for seg in range(4):
                        total += (
                            segregation_tensor[geno_sire, geno_dam, geno_child, seg]
                            * current_seg[seg, locus]
                        )
                    output[geno_sire, geno_dam, geno_child, locus] = total

    return output


@jit(nopython=True, nogil=True)
def project_child_genotypes(child_segs, child_values, output, n_loci):
    """Estimate the parental genotypes based on the child's genotypes and their segregation tensor.

    :param child_segs: the probability of each combination of
        the sire's genotype, the dam's genotype and the child's genotype of
        each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type child_segs: 4D numpy array of float32 with size 4 x 4 x 4 x n_loci
    :param child_values: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type child_values: 2D numpy array of float32 with size 4 x n_loci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype
        given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type output: 3D numpy array of float32 with size 4 x 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abci, ci -> abi", child_segs, child_values)
    for geno_sire in range(4):
        for geno_dam in range(4):
            for locus in range(n_loci):
                total = 0
                for geno_child in range(4):
                    total += (
                        child_segs[geno_sire, geno_dam, geno_child, locus]
                        * child_values[geno_child, locus]
                    )
                output[geno_sire, geno_dam, locus] = total

    return output


@jit(nopython=True, nogil=True)
def project_parent_genotypes(child_segs, parent_values, output, n_loci):
    """Project the parent genotypes down onto the child genotypes.

    :param child_segs: the probability of each combination of
        the sire's genotype, the dam's genotype and the child's genotype
        of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type child_segs: 4D numpy array of float32 with size 4 x 4 x 4 x n_loci
    :param parent_values: the probability of each combination of
        sire's genotype and the dams's genotype of each locus given the
        information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parent_values: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param output: the probability of child's genotype of each locus
        given the later and current generations from previous peeling cycle
        without the information of the current child (P(allele))
    :type output: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abci, abi -> ci", child_segs, parent_values)
    for geno_child in range(4):
        for locus in range(n_loci):
            total = 0
            for geno_sire in range(4):
                for geno_dam in range(4):
                    total += (
                        child_segs[geno_sire, geno_dam, geno_child, locus]
                        * parent_values[geno_sire, geno_dam, locus]
                    )
            output[geno_child, locus] = total

    return output


@jit(nopython=True, nogil=True)
def add_log_child_to_parents(child_to_parents, all_to_parents, n_loci):
    """Add one child's parent projection to the log-scale family total."""

    for geno_sire in range(4):
        for geno_dam in range(4):
            for locus in range(n_loci):
                all_to_parents[geno_sire, geno_dam, locus] += np.log(
                    child_to_parents[geno_sire, geno_dam, locus]
                )


@jit(nopython=True, nogil=True)
def add_log_child_to_parents_and_subtract_current(
    child_to_parents, all_to_parents, parents_minus_current_child, n_loci
):
    """Add one child to the family total and store its negative log contribution."""

    for geno_sire in range(4):
        for geno_dam in range(4):
            for locus in range(n_loci):
                log_value = np.log(child_to_parents[geno_sire, geno_dam, locus])
                all_to_parents[geno_sire, geno_dam, locus] += log_value
                parents_minus_current_child[geno_sire, geno_dam, locus] = -log_value


@jit(nopython=True, nogil=True)
def add_joint_parents_and_all_to_minus(
    parents_minus_child, joint_parents, all_to_parents, n_offspring, n_loci
):
    """Complete parent-minus-child log terms after all children are accumulated."""

    for child_index in range(n_offspring):
        for geno_sire in range(4):
            for geno_dam in range(4):
                for locus in range(n_loci):
                    parents_minus_child[child_index, geno_sire, geno_dam, locus] += (
                        np.log(joint_parents[geno_sire, geno_dam, locus])
                        + all_to_parents[geno_sire, geno_dam, locus]
                    )


@jit(nopython=True, nogil=True)
# pylint: disable=too-many-nested-blocks
def estimate_segregation_with_norm(
    segregation_tensors,
    genotype_values,
    output,
    n_loci,
):
    """Estimate segregation probabilities with tensor-specific normalization.

    :param segregation_tensor: the probability of each combination of
        the sire's genotype, the dam's genotype and the child's genotype
        and segregation without any other information (P(p, m, allele, seg))
    :type segregation_tensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param segregation_tensor_norm: the mean probability of each combination of
        the sire's genotype, the dam's genotype and the child's genotype across
        the child's segregation without any other information
        (1 / 4 x (P(p, m, allele)))
    :type segregation_tensor_norm: 3D numpy array of float32 with size 4 x 4 x 4
    :param parent_values: the probability of each combination of
        sire's genotype and the dams's genotype of each locus given the
        information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parent_values: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param child_values: the probability of each genotype of each locus of
        the current child given the information of itself and its later
        generations from previous peeling cycle (P(allele))
    :type child_values: 2D numpy array of float32 with size 4 x n_loci
    :param output: the probability of each segregation states of each
        locus of the current child given the information of later and
        current generations from previous peeling cycle (P(seg))
    :type output: 2D numpy array of float32 with size 4 x n_loci
    """
    segregation_tensor, segregation_tensor_norm = segregation_tensors
    parent_values, child_values = genotype_values
    # Equivalent einsum before normalization: output =
    # np.einsum("abcd, abi, ci -> di", segregation_tensor,
    # parent_values, child_values)
    for seg in range(4):
        for locus in range(n_loci):
            total = 0
            for geno_sire in range(4):
                for geno_dam in range(4):
                    for geno_child in range(4):
                        if (
                            segregation_tensor_norm[geno_sire, geno_dam, geno_child]
                            != 0
                        ):
                            total += (
                                segregation_tensor[geno_sire, geno_dam, geno_child, seg]
                                * parent_values[geno_sire, geno_dam, locus]
                                * child_values[geno_child, locus]
                                / segregation_tensor_norm[
                                    geno_sire, geno_dam, geno_child
                                ]
                            )
            output[seg, locus] = total
    return output


@jit(nopython=True, nogil=True)
def combine_and_reduce_axis1(joint_estimate, parent_estimate, output, n_loci):
    """Summing over axis 1 of joint_estimate with weights given by parent_estimate

    :param joint_estimate: the probability of each combination of
        sire's genotype and the dams's genotype of each locus
        given the information of later and current generations
        from previous peeling cycle (P(p, m))
    :type joint_estimate: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param parent_estimate: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type parent_estimate: 2D numpy array of float32 with size 4 x n_loci
    :return: the probability of each genotype of each locus of the sire
        given the information of later and current generations from previous peeling cycle (P(p))
    :rtype: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abi, bi -> ai", joint_estimate, parent_estimate)
    for geno_sire in range(4):
        for locus in range(n_loci):
            total = 0
            for geno_dam in range(4):
                total += (
                    joint_estimate[geno_sire, geno_dam, locus]
                    * parent_estimate[geno_dam, locus]
                )
            output[geno_sire, locus] = total
    return output


@jit(nopython=True, nogil=True)
def combine_and_reduce_axis0(joint_estimate, parent_estimate, output, n_loci):
    """Summing over axis 0 of joint_estimate with weights given by parent_estimate

    :param joint_estimate: the probability of each combination of
        sire's genotype and the dams's genotype of each locus
        given the information of later and current generations
        from previous peeling cycle (P(p, m))
    :type joint_estimate: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param parent_estimate: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type parent_estimate: 2D numpy array of float32 with size 4 x n_loci
    :return: the probability of each genotype of each locus of the dam
        given the information of later and current generations from previous peeling cycle (P(m))
    :rtype: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abi, ai -> bi", joint_estimate, parent_estimate)
    for geno_dam in range(4):
        for locus in range(n_loci):
            total = 0
            for geno_sire in range(4):
                total += (
                    joint_estimate[geno_sire, geno_dam, locus]
                    * parent_estimate[geno_sire, locus]
                )
            output[geno_dam, locus] = total
    return output


@jit(nopython=True, nogil=True)
def normalize_4_by_locus(mat, n_loci):
    """Normalize a 4 x n_loci probability matrix in place."""

    for locus in range(n_loci):
        total = 0
        for a in range(4):
            total += mat[a, locus]
        for a in range(4):
            mat[a, locus] /= total
    return mat


@jit(nopython=True, nogil=True)
def normalize_16_by_locus(mat, n_loci):
    """Normalize a 4 x 4 x n_loci probability tensor in place."""

    for locus in range(n_loci):
        total = 0
        for a in range(4):
            for b in range(4):
                total += mat[a, b, locus]
        for a in range(4):
            for b in range(4):
                mat[a, b, locus] /= total
    return mat


@jit(nopython=True, nogil=True)
def smooth_4_by_locus(mat, scale, floor, n_loci):
    """Apply a small in-place probability floor to a 4 x n_loci matrix."""

    for a in range(4):
        for locus in range(n_loci):
            mat[a, locus] = mat[a, locus] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def smooth_16_by_locus(mat, scale, floor, n_loci):
    """Apply a small in-place probability floor to a 4 x 4 x n_loci tensor."""

    for a in range(4):
        for b in range(4):
            for locus in range(n_loci):
                mat[a, b, locus] = mat[a, b, locus] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def exp_norm_2d(mat, n_loci):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 3D matrix with last axis represents the locus
    :type mat: 3D numpy array of float32 with size 4 x 4 x n_loci
    :return: the normalized exponential of the `mat`
    :rtype: 3D numpy array of float32 with size 4 x 4 x n_loci
    """
    for locus in range(n_loci):
        max_val = 1  # Log probabilities are non-positive, so 1 marks "unset".
        for a in range(4):
            for b in range(4):
                if mat[a, b, locus] > max_val or max_val == 1:
                    max_val = mat[a, b, locus]
        for a in range(4):
            for b in range(4):
                mat[a, b, locus] -= max_val
    tmp = np.exp(mat)
    normalize_16_by_locus(tmp, n_loci)
    return tmp


@jit(nopython=True, nogil=True)
def exp_norm_1d(mat, n_loci):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 2D matrix with last axis represents the locus
    :type mat: 2D numpy array of float32 with size 4 x n_loci
    :return: the normalized exponential of the `mat`
    :rtype: 2D numpy array of float32 with size 4 x n_loci
    """
    for locus in range(n_loci):
        max_val = 1  # Log probabilities are non-positive, so 1 marks "unset".
        for a in range(4):
            if mat[a, locus] > max_val or max_val == 1:
                max_val = mat[a, locus]
        for a in range(4):
            mat[a, locus] -= max_val
    tmp = np.exp(mat)
    normalize_4_by_locus(tmp, n_loci)
    return tmp


@jit(nopython=True, nogil=True)
def exp_norm_1d_in_place(mat, n_loci):
    """Take the exponential of a 4 x n_loci log matrix and normalize in place."""

    for locus in range(n_loci):
        max_val = 1
        for a in range(4):
            if mat[a, locus] > max_val or max_val == 1:
                max_val = mat[a, locus]
        total = 0
        for a in range(4):
            mat[a, locus] = np.exp(mat[a, locus] - max_val)
            total += mat[a, locus]
        for a in range(4):
            mat[a, locus] /= total
    return mat


@jit(nopython=True, nogil=True)
def collapse_segregation_in_place(segregation, transmission, forward, n_loci):
    """Using Baum-Welch algorithm to calculate the segregation probabilities.

    :param segregation: the probability of each segregation states of
        each locus of the current child given the information of later
        and current generations from previous peeling cycle
        (P(seg)). This is updated in place with the collapsed probabilities.
        the segregation state ordering: pp, pm, mp, mm
    :type segregation: 2D numpy array of float32 with size 4 x n_loci
    :param transmission: transmission function based on the distance between adjacent loci
    :type transmission: 1D numpy array of float32 with size (n_loci - 1)
    :param forward: work array for forward messages.
    :type forward: 2D numpy array of float32 with size 4 x n_loci
    :return: None. The function updates segregation in place with the collapsed probabilities.
    """
    # Forward-backward pass over segregation states ordered as pp, pm, mp, mm.
    tmp = np.full((4), 0, dtype=np.float32)
    new = np.full((4), 0, dtype=np.float32)
    prev = np.full((4), 1, dtype=np.float32)
    buffers = (tmp, new, prev)

    for state in range(4):
        forward[state, 0] = 1

    run_forward_segregation_pass(segregation, transmission, forward, buffers, n_loci)
    run_backward_segregation_pass(segregation, transmission, forward, buffers, n_loci)
    normalize_first_segregation_locus(segregation, forward, prev)


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e2": float32, "e1e": float32, "e2i": float32},
)
def run_forward_segregation_pass(segregation, transmission, forward, buffers, n_loci):
    """Run the forward pass for collapsed segregation probabilities."""

    tmp, new, prev = buffers
    for locus in range(1, n_loci):
        e = transmission[locus - 1]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2
        for state in range(4):
            tmp[state] = prev[state] * segregation[state, locus - 1]

        sum_state = 0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        for state in range(4):
            forward[state, locus] = new[state]
            prev[state] = new[state]


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e2": float32, "e1e": float32, "e2i": float32},
)
def run_backward_segregation_pass(segregation, transmission, forward, buffers, n_loci):
    """Run the backward pass for collapsed segregation probabilities."""

    tmp, new, prev = buffers
    for state in range(4):
        prev[state] = 1

    for locus in range(n_loci - 2, -1, -1):
        e = transmission[locus]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2

        for state in range(4):
            tmp[state] = prev[state] * segregation[state, locus + 1]

        sum_state = 0
        for state in range(4):
            sum_state += tmp[state]
        for state in range(4):
            tmp[state] = tmp[state] / sum_state

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        sum_state = 0
        for state in range(4):
            segregation[state, locus + 1] = (
                segregation[state, locus + 1] * forward[state, locus + 1] * prev[state]
            )
            sum_state += segregation[state, locus + 1]
        for state in range(4):
            segregation[state, locus + 1] = segregation[state, locus + 1] / sum_state
            prev[state] = new[state]


@jit(nopython=True, nogil=True)
def normalize_first_segregation_locus(segregation, forward, prev):
    """Normalize the first locus after the backward segregation pass."""

    sum_state = 0
    for state in range(4):
        segregation[state, 0] = segregation[state, 0] * forward[state, 0] * prev[state]
        sum_state += segregation[state, 0]
    for state in range(4):
        segregation[state, 0] = segregation[state, 0] / sum_state
