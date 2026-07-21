from numba import jit, float32
import numpy as np


# Numba handles integer state constants more reliably than string markers.
PEEL_UP = 0
PEEL_DOWN = 1


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e16": float32, "e1e": float32},
)
def peel_down(family, peeling_info, single_locus_mode):
    """Peel information down from parents to offspring.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`PeelingInfo.jit_peeling_information`
    :param single_locus_mode: A flag to indicate the mode of peeling.
        `False` if using multi-locus peeling, and `True` if using single-locus peeling.
    :type single_locus_mode: bool
    :return: None. The function modifies the peeling_info object in place.
    """
    is_x_chr = peeling_info.is_x_chr

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4
    e16 = e / 16

    anterior = peeling_info.anterior
    penetrance = peeling_info.penetrance
    posterior = peeling_info.posterior
    segregation = peeling_info.segregation

    segregation_tensor = peeling_info.segregation_tensor
    segregation_tensor_norm = peeling_info.segregation_tensor_norm

    n_loci = peeling_info.n_loci
    n_offspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    # Current child's projection onto joint parent genotypes.
    child_to_parents_current = np.full((4, 4, n_loci), 0, dtype=np.float32)
    # All children's projections accumulated on joint parent genotypes.
    all_to_parents = np.full((4, 4, n_loci), 0, dtype=np.float32)
    # Parent genotype probabilities after removing this family's current contribution.
    prob_sire = np.full((4, n_loci), 0, dtype=np.float32)
    prob_dam = np.full((4, n_loci), 0, dtype=np.float32)
    # Current child's posterior and penetrance terms collapsed by genotype state.
    child_values_current = np.full((4, n_loci), 0, dtype=np.float32)

    # Current child's segregation probabilities and forward-pass workspace.
    current_seg = np.full((4, n_loci), 1, dtype=np.float32)
    forward_seg = np.full((4, n_loci), 1, dtype=np.float32)

    needs_segregation_update = not single_locus_mode
    # Per-child inheritance tensors, reused later to update anterior probabilities.
    child_seg_tensor = np.full((n_offspring, 4, 4, 4, n_loci), 0, dtype=np.float32)
    # Joint parent estimates with one child's contribution removed.
    parents_minus_child = np.full((n_offspring, 4, 4, n_loci), 1, dtype=np.float32)
    if needs_segregation_update:
        # Stored child genotype values needed after parent-minus-child estimates are built.
        child_values_tensor = np.full((n_offspring, 4, n_loci), 0, dtype=np.float32)

    # Build parent genotype probabilities excluding this family's current contribution.
    setup_parent_genotype_probs(
        anterior,
        penetrance,
        posterior,
        peeling_info.posterior_sire_contribution[fam, :, :],
        peeling_info.posterior_dam_contribution[fam, :, :],
        sire,
        dam,
        prob_sire,
        prob_dam,
        n_loci,
    )

    # Equivalent einsum: joint_parents = np.einsum("ai, bi -> abi", prob_sire, prob_dam)
    joint_parents = get_joint_parents(prob_sire, prob_dam, n_loci)
    smooth_16_by_locus(joint_parents, e1e, e16, n_loci)

    if is_x_chr:
        for index in range(n_offspring):
            child = family.offspring[index]
            child_segs = child_seg_tensor[index, :, :, :, :]
            project_child_for_peel_x_chr(
                posterior,
                penetrance,
                segregation,
                peeling_info.sex,
                peeling_info.segregation_tensor_xy,
                peeling_info.segregation_tensor_xx,
                child,
                child_values_current,
                current_seg,
                child_segs,
                child_to_parents_current,
                n_loci,
                e1e,
                e4,
            )
            if needs_segregation_update:
                child_values_tensor[index, :, :] = child_values_current

            add_log_child_to_parents_and_subtract_current(
                child_to_parents_current,
                all_to_parents,
                parents_minus_child[index, :, :, :],
                n_loci,
            )
    else:
        for index in range(n_offspring):
            child = family.offspring[index]

            # Equivalent einsum: child_seg_tensor[index, :, :, :, :] = np.einsum("abcd, di -> abci", segregation_tensor, current_seg)
            child_segs = child_seg_tensor[index, :, :, :, :]

            # Equivalent einsum: child_to_parents_current = np.einsum("abci, ci -> abi", child_segs, child_values_current)
            project_child_for_peel_autosome(
                posterior,
                penetrance,
                segregation,
                peeling_info.segregation_tensor,
                child,
                child_values_current,
                current_seg,
                child_segs,
                child_to_parents_current,
                n_loci,
                e1e,
                e4,
            )
            if needs_segregation_update:
                child_values_tensor[index, :, :] = child_values_current

            add_log_child_to_parents_and_subtract_current(
                child_to_parents_current,
                all_to_parents,
                parents_minus_child[index, :, :, :],
                n_loci,
            )

    # Combine all child projections on the log scale, then subtract each child
    # to produce child-specific parent estimates.
    add_joint_parents_and_all_to_minus(
        parents_minus_child, joint_parents, all_to_parents, n_offspring, n_loci
    )

    # Convert from log scale to normalized probabilities.
    all_to_parents = exp_norm_2d(all_to_parents, n_loci)
    for i in range(n_offspring):
        parents_minus_child[i, :, :, :] = exp_norm_2d(
            parents_minus_child[i, :, :, :], n_loci
        )

    for i in range(n_offspring):
        child = family.offspring[i]

        # Equivalent einsum: anterior[child, :, :] = np.einsum("abci, abi -> ci", child_seg_tensor[i, :, :, :, :], parents_minus_child[i, :, :, :])
        child_anterior = anterior[child, :, :]
        project_parent_genotypes(
            child_seg_tensor[i, :, :, :, :],
            parents_minus_child[i, :, :, :],
            child_anterior,
            n_loci,
        )
        normalize_4_by_locus(child_anterior, n_loci)

    if needs_segregation_update:
        if is_x_chr:
            for i in range(n_offspring):
                child = family.offspring[i]
                child_values = child_values_tensor[i, :, :]

                if peeling_info.sex[child] == 0:  # 0=male, 1=female.
                    segregation_tensor = peeling_info.segregation_tensor_xy
                    segregation_tensor_norm = peeling_info.segregation_tensor_xy_norm
                elif peeling_info.sex[child] == 1:  # 0=male, 1=female.
                    segregation_tensor = peeling_info.segregation_tensor_xx
                    segregation_tensor_norm = peeling_info.segregation_tensor_xx_norm

                # Equivalent einsum: segregation[child, :, :] = np.einsum("abcd, abi, ci -> di", segregation_tensor, parents_minus_child[i, :, :, :], child_values)
                estimate_segregation_with_norm(
                    segregation_tensor,
                    segregation_tensor_norm,
                    parents_minus_child[i, :, :, :],
                    child_values,
                    segregation[child, :, :],
                    n_loci,
                )

                collapse_segregation_in_place(
                    segregation[child, :, :],
                    peeling_info.transmission_rate,
                    forward_seg,
                    n_loci,
                )
                for locus in range(n_loci):
                    for state in range(4):
                        segregation[child, state, locus] = (
                            e1e * segregation[child, state, locus] + e4
                        )
        else:
            for i in range(n_offspring):
                child = family.offspring[i]
                child_values = child_values_tensor[i, :, :]

                # Equivalent einsum: segregation[child, :, :] = np.einsum("abcd, abi, ci -> di", segregation_tensor, parents_minus_child[i, :, :, :], child_values)
                estimate_segregation_with_norm(
                    segregation_tensor,
                    segregation_tensor_norm,
                    parents_minus_child[i, :, :, :],
                    child_values,
                    segregation[child, :, :],
                    n_loci,
                )

                collapse_segregation_in_place(
                    segregation[child, :, :],
                    peeling_info.transmission_rate,
                    forward_seg,
                    n_loci,
                )
                for locus in range(n_loci):
                    for state in range(4):
                        segregation[child, state, locus] = (
                            e1e * segregation[child, state, locus] + e4
                        )


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e4": float32, "e1e": float32},
)
def peel_up(family, peeling_info):
    """Peel information up from offspring to parents.

    :param family: The family object that the peeling is performed on.
    :type family: class:`tinyhouse.Pedigree.Family`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`PeelingInfo.jit_peeling_information`
    :return: None. The function modifies the peeling_info object in place.
    """
    is_x_chr = peeling_info.is_x_chr

    e = 0.000001
    e1e = 1 - e
    e4 = e / 4

    anterior = peeling_info.anterior
    penetrance = peeling_info.penetrance
    posterior = peeling_info.posterior
    segregation = peeling_info.segregation

    n_loci = peeling_info.n_loci
    n_offspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    # Current child's projection onto joint parent genotypes.
    child_to_parents_current = np.full((4, 4, n_loci), 0, dtype=np.float32)
    # All children's projections accumulated on joint parent genotypes.
    all_to_parents = np.full((4, 4, n_loci), 0, dtype=np.float32)
    # Parent genotype probabilities after removing this family's current contribution.
    prob_sire = np.full((4, n_loci), 0, dtype=np.float32)
    prob_dam = np.full((4, n_loci), 0, dtype=np.float32)
    # Current child's posterior and penetrance terms collapsed by genotype state.
    child_values_current = np.full((4, n_loci), 0, dtype=np.float32)
    # Current child's segregation probabilities.
    current_seg = np.full((4, n_loci), 1, dtype=np.float32)
    # Single-child inheritance tensor reused for each offspring during peel-up.
    child_seg_tensor = np.full((1, 4, 4, 4, n_loci), 0, dtype=np.float32)
    child_segs_current = child_seg_tensor[0, :, :, :, :]

    setup_parent_genotype_probs(
        anterior,
        penetrance,
        posterior,
        peeling_info.posterior_sire_contribution[fam, :, :],
        peeling_info.posterior_dam_contribution[fam, :, :],
        sire,
        dam,
        prob_sire,
        prob_dam,
        n_loci,
    )
    smooth_4_by_locus(prob_sire, e1e, e4, n_loci)
    smooth_4_by_locus(prob_dam, e1e, e4, n_loci)

    if is_x_chr:
        for index in range(n_offspring):
            child = family.offspring[index]

            project_child_for_peel_x_chr(
                posterior,
                penetrance,
                segregation,
                peeling_info.sex,
                peeling_info.segregation_tensor_xy,
                peeling_info.segregation_tensor_xx,
                child,
                child_values_current,
                current_seg,
                child_segs_current,
                child_to_parents_current,
                n_loci,
                e1e,
                e4,
            )
            add_log_child_to_parents(child_to_parents_current, all_to_parents, n_loci)
    else:
        for index in range(n_offspring):
            child = family.offspring[index]

            project_child_for_peel_autosome(
                posterior,
                penetrance,
                segregation,
                peeling_info.segregation_tensor,
                child,
                child_values_current,
                current_seg,
                child_segs_current,
                child_to_parents_current,
                n_loci,
                e1e,
                e4,
            )
            add_log_child_to_parents(child_to_parents_current, all_to_parents, n_loci)

    all_to_parents = exp_norm_2d(all_to_parents, n_loci)

    sire_posterior = peeling_info.posterior_sire_contribution[fam, :, :]
    combine_and_reduce_axis1(all_to_parents, prob_dam, sire_posterior, n_loci)
    normalize_4_by_locus(sire_posterior, n_loci)
    smooth_4_by_locus(sire_posterior, e1e, e4, n_loci)

    dam_posterior = peeling_info.posterior_dam_contribution[fam, :, :]
    combine_and_reduce_axis0(all_to_parents, prob_sire, dam_posterior, n_loci)
    normalize_4_by_locus(dam_posterior, n_loci)
    smooth_4_by_locus(dam_posterior, e1e, e4, n_loci)


# JIT helpers below replace the einsum calls used in the original implementation.


@jit(nopython=True, nogil=True)
def setup_parent_genotype_probs(
    anterior,
    penetrance,
    posterior,
    sire_contribution,
    dam_contribution,
    sire,
    dam,
    prob_sire,
    prob_dam,
    n_loci,
):
    """Build the normalized sire and dam genotype probabilities for a family."""

    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] = np.log(posterior[sire, state, locus]) - np.log(
                sire_contribution[state, locus]
            )
            prob_dam[state, locus] = np.log(posterior[dam, state, locus]) - np.log(
                dam_contribution[state, locus]
            )

    exp_norm_1d_in_place(prob_sire, n_loci)
    exp_norm_1d_in_place(prob_dam, n_loci)

    for state in range(4):
        for locus in range(n_loci):
            prob_sire[state, locus] *= (
                anterior[sire, state, locus] * penetrance[sire, state, locus]
            )
            prob_dam[state, locus] *= (
                anterior[dam, state, locus] * penetrance[dam, state, locus]
            )

    normalize_4_by_locus(prob_sire, n_loci)
    normalize_4_by_locus(prob_dam, n_loci)


@jit(
    nopython=True,
    nogil=True,
    locals={"e4": float32, "e1e": float32},
)
def project_child_for_peel_autosome(
    posterior,
    penetrance,
    segregation,
    segregation_tensor,
    child,
    child_values,
    current_seg,
    child_segs,
    child_to_parents,
    n_loci,
    e1e,
    e4,
):
    """Project one autosomal child onto parental genotypes."""

    for state in range(4):
        for locus in range(n_loci):
            child_values[state, locus] = (
                posterior[child, state, locus] * penetrance[child, state, locus]
            )
            current_seg[state, locus] = segregation[child, state, locus]

    normalize_4_by_locus(child_values, n_loci)
    smooth_4_by_locus(child_values, e1e, e4, n_loci)
    normalize_4_by_locus(current_seg, n_loci)

    create_child_segs(segregation_tensor, current_seg, child_segs, n_loci)
    project_child_genotypes(child_segs, child_values, child_to_parents, n_loci)


@jit(
    nopython=True,
    nogil=True,
    locals={"e4": float32, "e1e": float32},
)
def project_child_for_peel_x_chr(
    posterior,
    penetrance,
    segregation,
    sex,
    segregation_tensor_xy,
    segregation_tensor_xx,
    child,
    child_values,
    current_seg,
    child_segs,
    child_to_parents,
    n_loci,
    e1e,
    e4,
):
    """Project one X-chromosome child onto parental genotypes."""

    for state in range(4):
        for locus in range(n_loci):
            child_values[state, locus] = (
                posterior[child, state, locus] * penetrance[child, state, locus]
            )
            current_seg[state, locus] = segregation[child, state, locus]

    normalize_4_by_locus(child_values, n_loci)
    smooth_4_by_locus(child_values, e1e, e4, n_loci)
    normalize_4_by_locus(current_seg, n_loci)

    if sex[child] == 0:  # 0=male, 1=female.
        segregation_tensor = segregation_tensor_xy
    else:
        segregation_tensor = segregation_tensor_xx

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
    for a in range(4):
        for b in range(4):
            for i in range(n_loci):
                output[a, b, i] = prob_sire[a, i] * prob_dam[b, i]
    return output


@jit(nopython=True, nogil=True)
def create_child_segs(segregation_tensor, current_seg, output, n_loci):
    """Creates the child-specific segregation tensor using the child's current segregation estimate.

    :param segregation_tensor: the probability of each combination of the sire's genotype, the dam's genotype,
        child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregation_tensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param current_seg: The probability of each segregation of each locus of the child
        with information from previous peeling cycle (P(seg))
    :type current_seg: 2D numpy array of float32 with size 4 x n_loci
    :param output: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle (P(p, m, allele))
    :type output: 4D numpy array of float32 with size 4 x 4 x 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abcd, di -> abci", segregation_tensor, current_seg)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for i in range(n_loci):
                    total = 0
                    for d in range(4):
                        total += segregation_tensor[a, b, c, d] * current_seg[d, i]
                    output[a, b, c, i] = total

    return output


@jit(nopython=True, nogil=True)
def project_child_genotypes(child_segs, child_values, output, n_loci):
    """Estimate the parental genotypes based on the child's genotypes and their segregation tensor.

    :param child_segs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
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
    for a in range(4):
        for b in range(4):
            for i in range(n_loci):
                total = 0
                for c in range(4):
                    total += child_segs[a, b, c, i] * child_values[c, i]
                output[a, b, i] = total
    return output


@jit(nopython=True, nogil=True)
def project_parent_genotypes(child_segs, parent_values, output, n_loci):
    """Project the parent genotypes down onto the child genotypes.

    :param child_segs: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype of each locus with information from previous peeling cycle
        (P(p, m, allele))
    :type child_segs: 4D numpy array of float32 with size 4 x 4 x 4 x n_loci
    :param parent_values: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parent_values: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param output: the probability of child's genotype of each locus given the later and current generations
        from previous peeling cycle without the information of the current child (P(allele))
    :type output: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abci, abi -> ci", child_segs, parent_values)

    for c in range(4):
        for i in range(n_loci):
            total = 0
            for a in range(4):
                for b in range(4):
                    total += child_segs[a, b, c, i] * parent_values[a, b, i]
            output[c, i] = total

    return output


@jit(nopython=True, nogil=True)
def add_log_child_to_parents(child_to_parents, all_to_parents, n_loci):
    """Add one child's parent projection to the log-scale family total."""

    for a in range(4):
        for b in range(4):
            for i in range(n_loci):
                all_to_parents[a, b, i] += np.log(child_to_parents[a, b, i])


@jit(nopython=True, nogil=True)
def add_log_child_to_parents_and_subtract_current(
    child_to_parents, all_to_parents, parents_minus_current_child, n_loci
):
    """Add one child to the family total and store its negative log contribution."""

    for a in range(4):
        for b in range(4):
            for i in range(n_loci):
                log_value = np.log(child_to_parents[a, b, i])
                all_to_parents[a, b, i] += log_value
                parents_minus_current_child[a, b, i] = -log_value


@jit(nopython=True, nogil=True)
def add_joint_parents_and_all_to_minus(
    parents_minus_child, joint_parents, all_to_parents, n_offspring, n_loci
):
    """Complete parent-minus-child log terms after all children are accumulated."""

    for child_index in range(n_offspring):
        for a in range(4):
            for b in range(4):
                for i in range(n_loci):
                    parents_minus_child[child_index, a, b, i] += (
                        np.log(joint_parents[a, b, i]) + all_to_parents[a, b, i]
                    )


@jit(nopython=True, nogil=True)
def estimate_segregation_with_norm(
    segregation_tensor,
    segregation_tensor_norm,
    parent_values,
    child_values,
    output,
    n_loci,
):
    """Estimate segregation probabilities with tensor-specific normalization.

    :param segregation_tensor: the probability of each combination of the sire's genotype, the dam's genotype and
        the child's genotype and segregation without any other information (P(p, m, allele, seg))
    :type segregation_tensor: 4D numpy array of float32 with size 4 x 4 x 4 x 4
    :param segregation_tensor_norm: the mean probability of each combination of the sire's genotype, the dam's genotype
        and the child's genotype across the child's segregation without any other information
        (1 / 4 x (P(p, m, allele)))
    :type segregation_tensor_norm: 3D numpy array of float32 with size 4 x 4 x 4
    :param parent_values: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        without the information of the current child (P(p, m))
    :type parent_values: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param child_values: the probability of each genotype of each locus of the current child
        given the information of itself and its later generations from previous peeling cycle
        (P(allele))
    :type child_values: 2D numpy array of float32 with size 4 x n_loci
    :param output: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
        (P(seg))
    :type output: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum before normalization: output = np.einsum("abcd, abi, ci -> di", segregation_tensor, parent_values, child_values)
    for d in range(4):
        for i in range(n_loci):
            total = 0
            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        if segregation_tensor_norm[a, b, c] != 0:
                            total += (
                                segregation_tensor[a, b, c, d]
                                * parent_values[a, b, i]
                                * child_values[c, i]
                                / segregation_tensor_norm[a, b, c]
                            )
            output[d, i] = total
    return output


@jit(nopython=True, nogil=True)
def combine_and_reduce_axis1(joint_estimate, parent_estimate, output, n_loci):
    """Summing over axis 1 of joint_estimate with weights given by parent_estimate

    :param joint_estimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type joint_estimate: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param parent_estimate: the probability of each genotype of each locus of the dam
        with information of the current child from previous peeling cycle (P(m))
    :type parent_estimate: 2D numpy array of float32 with size 4 x n_loci
    :return: the probability of each genotype of each locus of the sire
        given the information of later and current generations from previous peeling cycle (P(p))
    :rtype: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abi, bi -> ai", joint_estimate, parent_estimate)
    for a in range(4):
        for i in range(n_loci):
            total = 0
            for b in range(4):
                total += joint_estimate[a, b, i] * parent_estimate[b, i]
            output[a, i] = total
    return output


@jit(nopython=True, nogil=True)
def combine_and_reduce_axis0(joint_estimate, parent_estimate, output, n_loci):
    """Summing over axis 0 of joint_estimate with weights given by parent_estimate

    :param joint_estimate: the probability of each combination of sire's genotype and the dams's genotype
        of each locus given the information of later and current generations from previous peeling cycle
        (P(p, m))
    :type joint_estimate: 3D numpy array of float32 with size 4 x 4 x n_loci
    :param parent_estimate: the probability of each genotype of each locus of the sire
        with information of the current child from previous peeling cycle (P(p))
    :type parent_estimate: 2D numpy array of float32 with size 4 x n_loci
    :return: the probability of each genotype of each locus of the dam
        given the information of later and current generations from previous peeling cycle (P(m))
    :rtype: 2D numpy array of float32 with size 4 x n_loci
    """
    # Equivalent einsum: output = np.einsum("abi, ai -> bi", joint_estimate, parent_estimate)
    for b in range(4):
        for i in range(n_loci):
            total = 0
            for a in range(4):
                total += joint_estimate[a, b, i] * parent_estimate[a, i]
            output[b, i] = total
    return output


@jit(nopython=True, nogil=True)
def normalize_4_by_locus(mat, n_loci):
    """Normalize a 4 x n_loci probability matrix in place."""

    for i in range(n_loci):
        total = 0
        for a in range(4):
            total += mat[a, i]
        for a in range(4):
            mat[a, i] /= total
    return mat


@jit(nopython=True, nogil=True)
def normalize_16_by_locus(mat, n_loci):
    """Normalize a 4 x 4 x n_loci probability tensor in place."""

    for i in range(n_loci):
        total = 0
        for a in range(4):
            for b in range(4):
                total += mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] /= total
    return mat


@jit(nopython=True, nogil=True)
def smooth_4_by_locus(mat, scale, floor, n_loci):
    """Apply a small in-place probability floor to a 4 x n_loci matrix."""

    for a in range(4):
        for i in range(n_loci):
            mat[a, i] = mat[a, i] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def smooth_16_by_locus(mat, scale, floor, n_loci):
    """Apply a small in-place probability floor to a 4 x 4 x n_loci tensor."""

    for a in range(4):
        for b in range(4):
            for i in range(n_loci):
                mat[a, b, i] = mat[a, b, i] * scale + floor
    return mat


@jit(nopython=True, nogil=True)
def exp_norm_2d(mat, n_loci):
    """Output is to take the exponential of the matrix and normalize each locus.

    :param mat: a 3D matrix with last axis represents the locus
    :type mat: 3D numpy array of float32 with size 4 x 4 x n_loci
    :return: the normalized exponential of the `mat`
    :rtype: 3D numpy array of float32 with size 4 x 4 x n_loci
    """
    for i in range(n_loci):
        max_val = 1  # Log probabilities are non-positive, so 1 marks "unset".
        for a in range(4):
            for b in range(4):
                if mat[a, b, i] > max_val or max_val == 1:
                    max_val = mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] -= max_val
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
    for i in range(n_loci):
        max_val = 1  # Log probabilities are non-positive, so 1 marks "unset".
        for a in range(4):
            if mat[a, i] > max_val or max_val == 1:
                max_val = mat[a, i]
        for a in range(4):
            mat[a, i] -= max_val
    tmp = np.exp(mat)
    normalize_4_by_locus(tmp, n_loci)
    return tmp


@jit(nopython=True, nogil=True)
def exp_norm_1d_in_place(mat, n_loci):
    """Take the exponential of a 4 x n_loci log matrix and normalize in place."""

    for i in range(n_loci):
        max_val = 1
        for a in range(4):
            if mat[a, i] > max_val or max_val == 1:
                max_val = mat[a, i]
        total = 0
        for a in range(4):
            mat[a, i] = np.exp(mat[a, i] - max_val)
            total += mat[a, i]
        for a in range(4):
            mat[a, i] /= total
    return mat


@jit(
    nopython=True,
    nogil=True,
    locals={"e": float32, "e2": float32, "e1e": float32, "e2i": float32},
)
def collapse_segregation_in_place(segregation, transmission, forward, n_loci):
    """Using Baum-Welch algorithm to calculate the segregation probabilities.

    :param segregation: the probability of each segregation states of each locus of the current child
        given the information of later and current generations from previous peeling cycle
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

    for j in range(4):
        forward[j, 0] = 1

    for i in range(1, n_loci):
        e = transmission[i - 1]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2
        for j in range(4):
            tmp[j] = prev[j] * segregation[j, i - 1]

        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j] / sum_j

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        for j in range(4):
            forward[j, i] = new[j]
            prev[j] = new[j]

    for j in range(4):
        prev[j] = 1

    for i in range(n_loci - 2, -1, -1):
        e = transmission[i]
        e2 = e**2
        e1e = e * (1 - e)
        e2i = (1.0 - e) ** 2

        for j in range(4):
            tmp[j] = prev[j] * segregation[j, i + 1]

        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j] / sum_j

        new[0] = e2 * tmp[3] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[0]
        new[1] = e2 * tmp[2] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[1]
        new[2] = e2 * tmp[1] + e1e * (tmp[0] + tmp[3]) + e2i * tmp[2]
        new[3] = e2 * tmp[0] + e1e * (tmp[1] + tmp[2]) + e2i * tmp[3]

        sum_j = 0
        for j in range(4):
            segregation[j, i + 1] = segregation[j, i + 1] * forward[j, i + 1] * prev[j]
            sum_j += segregation[j, i + 1]
        for j in range(4):
            segregation[j, i + 1] = segregation[j, i + 1] / sum_j
            prev[j] = new[j]

    sum_j = 0
    for j in range(4):
        segregation[j, 0] = segregation[j, 0] * forward[j, 0] * prev[j]
        sum_j += segregation[j, 0]
    for j in range(4):
        segregation[j, 0] = segregation[j, 0] / sum_j

    return
