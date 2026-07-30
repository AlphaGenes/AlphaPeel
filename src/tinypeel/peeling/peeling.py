"""Main peeling functions."""

import concurrent.futures

from numba import jit, float32
import numpy as np

from .peeling_kernels import (
    add_joint_parents_and_all_to_minus,
    add_log_child_to_parents,
    add_log_child_to_parents_and_subtract_current,
    collapse_segregation_in_place,
    combine_and_reduce_axis0,
    combine_and_reduce_axis1,
    create_child_segs,
    estimate_segregation_with_norm,
    exp_norm_1d,
    exp_norm_1d_in_place,
    exp_norm_2d,
    get_joint_parents,
    normalize_4_by_locus,
    project_child_genotypes,
    project_parent_genotypes,
    smooth_4_by_locus,
    smooth_16_by_locus,
)


def update_posterior(peeling_info, sires, dams):
    """Updates the posterior term for a specific set of sires and dams.

    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param sires: collection of sires to update
    :type sires: set of class:`tinyhouse.Pedigree.Individual`
    :param dams: collection of dams to update
    :type dams: set of class:`tinyhouse.Pedigree.Individual`
    :return: None. The function modifies the peeling_info object in place
    """

    for sire in sires:
        update_sire(sire, peeling_info)

    for dam in dams:
        update_dam(dam, peeling_info)


def update_sire(sire, peeling_info):
    """Updates the posterior term for a specific sire.

    :param sire: the sire to update
    :type sire: class: `tinyhouse.Pedigree.Individual`
    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function modifies the peeling_info object in place
    """
    fam_list = [fam.idn for fam in sire.families]
    sire = sire.idn
    sire_posterior = peeling_info.posterior[sire, :, :]
    sire_posterior[:, :] = 0
    for fam_id in fam_list:
        log_update = np.log(peeling_info.posterior_sire_contribution[fam_id, :, :])
        sire_posterior += log_update

    # Convert accumulated log terms back to normalized probabilities.
    sire_posterior[:, :] = exp_norm_1d(sire_posterior, peeling_info.n_loci)
    sire_posterior /= np.sum(sire_posterior, 0)


def update_dam(dam, peeling_info):
    """Updates the posterior term for a specific dam.

    :param dam: the dam to update
    :type dam: class: `tinyhouse.Pedigree.Individual`
    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function modifies the peeling_info object in place
    """
    fam_list = [fam.idn for fam in dam.families]
    dam = dam.idn
    dam_posterior = peeling_info.posterior[dam, :, :]
    dam_posterior[:, :] = 0
    for fam_id in fam_list:
        log_update = np.log(peeling_info.posterior_dam_contribution[fam_id, :, :])
        dam_posterior += log_update

    dam_posterior[:, :] = exp_norm_1d(dam_posterior, peeling_info.n_loci)
    dam_posterior /= np.sum(dam_posterior, 0)


def peel_down(family, peeling_info, single_locus_mode, locus_thread_blocks=None):
    """Peel information down from parents to offspring."""

    if locus_thread_blocks is None:
        _peel_down_serial(family, peeling_info, single_locus_mode)
        if not single_locus_mode:
            collapse_child_segregations_after_peel_down(
                family, peeling_info, 0.999999, 0.000001 / 4
            )
        return

    _peel_down_threaded_parent_setup(
        family, peeling_info, single_locus_mode, locus_thread_blocks
    )


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e16": float32, "e1e": float32},
)
def _peel_down_serial(family, peeling_info, single_locus_mode):
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
        estimate_child_segregations_for_peel_down(family, peeling_info, workspace)


def _peel_down_threaded_parent_setup(
    family, peeling_info, single_locus_mode, locus_thread_blocks
):
    """Peel down with each Python thread owning a view-backed locus block."""

    n_workers = len(locus_thread_blocks)
    futures = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for _, _, block_info in locus_thread_blocks:
            futures.append(
                executor.submit(
                    _peel_down_locus_block,
                    family,
                    block_info,
                    single_locus_mode,
                )
            )

        for future in futures:
            future.result()

    if not single_locus_mode:
        collapse_child_segregations_after_peel_down(
            family, peeling_info, 0.999999, 0.000001 / 4
        )


def _peel_down_locus_block(family, block_info, single_locus_mode):
    """Peel down one view-backed locus block."""

    _peel_down_serial(family, block_info, single_locus_mode)


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
def estimate_child_segregations_for_peel_down(family, peeling_info, workspace):
    """Estimate uncollapsed child segregation probabilities after peel-down."""

    for index, child in enumerate(family.offspring):
        estimate_one_child_segregation_for_peel_down(
            child, peeling_info, workspace, index
        )


@jit(nopython=True, nogil=True)
def collapse_child_segregations_after_peel_down(family, peeling_info, scale, floor):
    """Collapse and smooth child segregation probabilities after peel-down."""

    for child in family.offspring:
        collapse_one_child_segregation_for_peel_down(child, peeling_info)
        smooth_4_by_locus(
            peeling_info.segregation[child, :, :],
            scale,
            floor,
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
def collapse_one_child_segregation_for_peel_down(child, peeling_info):
    """Collapse one child's full-chromosome segregation probabilities."""

    forward = np.full((4, peeling_info.n_loci), 1, dtype=np.float32)
    collapse_segregation_in_place(
        peeling_info.segregation[child, :, :],
        peeling_info.transmission_rate,
        forward,
        peeling_info.n_loci,
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


def peel_up(family, peeling_info, locus_thread_blocks=None):
    """Peel information up from offspring to parents."""

    if locus_thread_blocks is None:
        _peel_up_serial(family, peeling_info)
        return

    _peel_up_threaded_by_locus(family, locus_thread_blocks)


def _peel_up_threaded_by_locus(family, locus_thread_blocks):
    """Peel up with each Python thread owning a view-backed locus block."""

    n_workers = len(locus_thread_blocks)
    futures = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
        for _, _, block_info in locus_thread_blocks:
            futures.append(
                executor.submit(
                    _peel_up_serial,
                    family,
                    block_info,
                )
            )

        for future in futures:
            future.result()


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
