"""Update allele frequencies, error rates, and phenotype penetrance estimates."""

import warnings

from numba import jit
import numpy as np

from .peeling_info_module import get_base_individual_penetrance, get_het_midpoint
from ..tinyhouse.ProbMath import (
    getGenotypesFromMultiMaf,
    getGenotypesFromMaf,
    updateGenoProbsFromPhenotype,
)


# Alternative allele frequencies are estimated with Newton updates under a
# Hardy-Weinberg genotype prior.


def prepare_alternative_allele_probabilities(pedigree, peeling_info, args):
    """Prepare initial alternative allele probabilities before peeling cycles."""

    if args.alt_allele_prob_file is not None:
        prepare_input_alternative_allele_probabilities(pedigree, peeling_info)
    else:
        prepare_default_alternative_allele_probabilities(pedigree, peeling_info)

    if args.est_start_alt_allele_prob:
        if args.alt_allele_prob_file is not None and len(pedigree.AAP) > 1:
            warnings.warn(
                "-est_start_alt_allele_prob will overwrite any differences "
                "between metafounders. "
                "To avoid this, please use -est_alt_allele_prob instead."
            )
        update_maf(pedigree, peeling_info)


def prepare_input_alternative_allele_probabilities(pedigree, peeling_info):
    """Prepare user-supplied alternative allele probabilities for used metafounders."""

    mf_pedigree = []
    maf_geno_cache = {}
    for ind in pedigree:
        if ind.isFounder() and ind.MetaFounder is not None:
            update_founder_metafounder_priors(ind, pedigree, peeling_info, mf_pedigree)
            maf_geno = get_maf_genotypes_for_meta_founder(
                ind.MetaFounder, pedigree, peeling_info.n_loci, maf_geno_cache
            )
            peeling_info.anterior[ind.idn, :, :] = maf_geno
    remove_unused_metafounder_priors(pedigree, mf_pedigree)


def update_founder_metafounder_priors(ind, pedigree, peeling_info, mf_pedigree):
    """Add and validate alternative allele priors for one founder's metafounders."""

    for mfx in ind.MetaFounder:
        if mfx in mf_pedigree:
            continue

        mf_pedigree.append(mfx)
        if pedigree.AAP.get(mfx) is None:
            pedigree.AAP[mfx] = np.full(peeling_info.n_loci, 0.5, dtype=np.float32)
        else:
            validate_and_clip_alternative_allele_probabilities(
                pedigree.AAP[mfx], mfx, peeling_info.n_loci
            )


def validate_and_clip_alternative_allele_probabilities(aap, mfx, n_loci):
    """Validate and clip one metafounder's alternative allele probabilities."""

    for i in range(n_loci):
        aap_value = aap[i]
        if aap_value > 1 or aap_value < 0:
            raise ValueError(
                f"Invalid value {aap_value} for alternative allele probability "
                f"for metafounder {mfx} at locus {i}. \n"
                "Values must be between 0 and 1. Set to 0.5 (default) if unknown."
            )
        if aap_value < 0.001:
            aap[i] = 0.001
        elif aap_value > 0.999:
            aap[i] = 0.999


def remove_unused_metafounder_priors(pedigree, mf_pedigree):
    """Drop alternative allele priors for metafounders absent from the pedigree."""

    mf_input = pedigree.AAP.copy()
    for mfx in mf_input:
        if mfx not in mf_pedigree:
            del pedigree.AAP[mfx]
            warnings.warn(
                f"{mfx} is not in the pedigree. "
                f"The alternative allele probability for {mfx} has been ignored."
            )


def prepare_default_alternative_allele_probabilities(pedigree, peeling_info):
    """Create default alternative allele probabilities for all metafounders."""

    for ind in pedigree:
        if ind.MetaFounder is not None:
            for mfx in ind.MetaFounder:
                if pedigree.AAP.get(mfx) is None:
                    pedigree.AAP[mfx] = np.full(
                        peeling_info.n_loci, 0.5, dtype=np.float32
                    )


def update_maf(pedigree, peeling_info):
    """Estimates the alternative allele frequency at all loci (i.e markers).

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function updates the pedigree.AAP attribute with
    the new alternative allele frequencies.
    """
    if peeling_info.is_x_chr:
        warnings.warn(
            "Updating error rates and alternative allele frequencies \
                for X chromosomes is not well tested. \
                Recommend running without that option."
        )
    meta_founders = list(pedigree.AAP.keys())
    genotyped_by_locus = get_genotyped_status(
        pedigree, peeling_info.n_ind, peeling_info.n_loci
    )
    for i in range(peeling_info.n_loci):
        genotyped = genotyped_by_locus[:, i]
        for mfx in meta_founders:
            alternative_allele_prob = pedigree.AAP[mfx]
            alternative_allele_prob[i] = newton_maf_updates(
                peeling_info, alternative_allele_prob, i, genotyped
            )

    for mfx in meta_founders:
        pedigree.AAP[mfx] = pedigree.AAP[mfx].astype(np.float32)

    maf_geno_cache = {}
    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            maf_geno = get_maf_genotypes_for_meta_founder(
                ind.MetaFounder, pedigree, peeling_info.n_loci, maf_geno_cache
            )
            peeling_info.anterior[ind.idn, :, :] = maf_geno


def individual_has_observed_data_at_locus(ind, index):
    """Return whether an individual has genotype or read data at one locus."""

    if ind.genotypes is not None and ind.genotypes[index] != 9:
        return True

    if ind.reads is not None:
        return ind.reads[0][index] != 0 or ind.reads[1][index] != 0

    return False


def get_genotyped_status_for_locus(pedigree, n_ind, index):
    """Build the observed-data status vector for one locus."""

    genotyped = np.full(n_ind, False, dtype=np.bool_)
    for ind in pedigree:
        genotyped[ind.idn] = individual_has_observed_data_at_locus(ind, index)

    return genotyped


def get_genotyped_status(pedigree, n_ind, n_loci):
    """Build the observed-data status matrix for all loci."""

    genotyped = np.full((n_ind, n_loci), False, dtype=np.bool_)
    for ind in pedigree:
        ind_genotyped = genotyped[ind.idn, :]
        if ind.genotypes is not None:
            ind_genotyped |= ind.genotypes != 9
        if ind.reads is not None:
            ind_genotyped |= (ind.reads[0] != 0) | (ind.reads[1] != 0)

    return genotyped


def get_maf_genotypes_for_meta_founder(meta_founder, pedigree, n_loci, cache):
    """Return the MAF genotype prior for a metafounder tuple."""

    key = tuple(meta_founder)
    if key not in cache:
        if len(meta_founder) == 2:
            alternative_allele_prob = {
                k: np.zeros(n_loci, dtype=np.float32) for k in meta_founder
            }
            for mfx in meta_founder:
                alternative_allele_prob[mfx] = pedigree.AAP[mfx]
            cache[key] = getGenotypesFromMultiMaf(alternative_allele_prob)
        else:
            cache[key] = getGenotypesFromMaf(pedigree.AAP[meta_founder[0]])

    return cache[key]


def newton_maf_updates(peeling_info, alternative_allele_prob, index, genotyped):
    """Iterative approximation for the prior alternative allele frequency.
    Currently limits all alternative_allele_prob to be between 0.001 and 0.999.

    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param alternative_allele_prob: starting alternative allele frequencies, default 0.5
    :type alternative_allele_prob: 1D numpy array with length equal to the number of loci
    :param index: the marker index for which to update the alternative allele frequency
    :type index: int
    :param genotyped: whether each individual has observed data for this locus
    :type genotyped: 1D numpy array of bool
    :return: the updated alternative allele frequency for the given marker index
    :rtype: float
    """

    if alternative_allele_prob[index] < 0.001:
        maf = 0.001
    elif alternative_allele_prob[index] > 0.999:
        maf = 0.999
    else:
        maf = alternative_allele_prob[index]

    iters = 5
    converged = False
    while not converged:
        maf_old = maf
        delta = get_newton_update(maf_old, peeling_info, index, genotyped)
        maf = maf_old + delta
        maf = max(maf, 0.001)
        maf = min(maf, 0.999)
        if abs(maf - maf_old) < 0.0001:
            converged = True
        iters -= 1
        if iters < 0:
            converged = True
    return maf


@jit(nopython=True)
def get_newton_update(p, peeling_info, index, genotyped):
    """Calculates the alternative allele frequency using Newton's method of optimisation.

    :param p: the current alternative allele frequency estimate
    :type p: float
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param index: the marker index for which to update the alternative allele frequency
    :type index: int
    :param genotyped: whether each individual has observed data for this locus
    :type genotyped: 1D numpy array of bool
    :return: ratio of the first and second derivatives of the log likelihood function
    to be added to the current alternative allele frequency estimate.
    :rtype: float
    """
    # First and second derivatives of the log likelihood.
    ll_p = 0
    ll_pp = 0

    # Add one pseudo-observation for each homozygous/heterozygous state.
    delta_p, delta_pp = add_individual_scalars_to_update(1, 0, 0, p)
    ll_p, ll_pp = add_update(ll_p, ll_pp, delta_p, delta_pp)
    delta_p, delta_pp = add_individual_scalars_to_update(0, 1, 0, p)
    ll_p, ll_pp = add_update(ll_p, ll_pp, 2 * delta_p, 2 * delta_pp)
    delta_p, delta_pp = add_individual_scalars_to_update(0, 0, 1, p)
    ll_p, ll_pp = add_update(ll_p, ll_pp, delta_p, delta_pp)
    for i in range(peeling_info.n_ind):
        if genotyped[i]:
            d0 = peeling_info.penetrance[i, 0, index]
            d1 = (
                peeling_info.penetrance[i, 1, index]
                + peeling_info.penetrance[i, 2, index]
            )
            d2 = peeling_info.penetrance[i, 3, index]
            delta_p, delta_pp = add_individual_scalars_to_update(d0, d1, d2, p)
            ll_p, ll_pp = add_update(ll_p, ll_pp, delta_p, delta_pp)
    # No observed data can leave the Newton derivative terms at zero.
    if ll_p == 0 or ll_pp == 0:
        return 0
    return -ll_p / ll_pp


@jit(nopython=True)
def add_update(ll_p, ll_pp, delta_p, delta_pp):
    """Adds the first and second derivatives of the log likelihood for one individual to the total.

    :param ll_p: the first derivative of the log likelihood
    :type ll_p: float
    :param ll_pp: the second derivative of the log likelihood
    :type ll_pp: float
    :param delta_p: the first derivative of the log likelihood for one individual
    :type delta_p: float
    :param delta_pp: the second derivative of the log likelihood for one individual
    :type delta_pp: float
    :return: updated first and second derivatives of the log likelihood
    :rtype: tuple(float, float)
    """
    ll_p += delta_p
    ll_pp += delta_pp
    return ll_p, ll_pp


@jit(nopython=True)
def add_individual_to_update(d, p, ll_p, ll_pp):
    """Adds each available genotype data to first and second derivatives of the log likelihood.

    :param d: the penetrance term for genotyped individuals as genotype probabilities
    :type d: 1D numpy array with length 4 (i.e [p(AA), p(aA), p(Aa), p(aa)])
    :param p: the current alternative allele frequency estimate
    :type p: float
    :param ll_p: the first derivative of the log likelihood
    :type ll_p: float
    :param ll_pp: the second derivative of the log likelihood
    :type ll_pp: float
    :return: updated first and second derivatives of the log likelihood
    :rtype: tuple(float, float)
    """
    d0 = d[0]
    d1 = d[1] + d[2]
    d2 = d[3]

    delta_p, delta_pp = add_individual_scalars_to_update(d0, d1, d2, p)
    ll_p, ll_pp = add_update(ll_p, ll_pp, delta_p, delta_pp)
    return ll_p, ll_pp


@jit(nopython=True)
def add_individual_scalars_to_update(d0, d1, d2, p):
    """Adds pre-collapsed genotype probabilities to the Newton update terms."""

    f = d0 * (1 - p) ** 2 + d1 * p * (1 - p) + d2 * p**2
    fp = (d1 - 2 * d0) + 2 * p * (d0 + d2 - d1)
    fpp = 2 * (d0 + d2 - d1)

    return fp / f, fpp / f - (fp / f) ** 2


def add_founder_to_alt_allele_update(ind, geno_probs, alternative_allele_prob, ind_mf):
    """Accumulate one founder's allele-probability contribution."""

    if len(ind.MetaFounder) == 2:
        add_two_metafounder_contribution(
            ind.MetaFounder,
            geno_probs,
            alternative_allele_prob,
            ind_mf,
        )
    else:
        add_single_metafounder_contribution(
            ind.MetaFounder[0],
            geno_probs,
            alternative_allele_prob,
            ind_mf,
        )


def add_two_metafounder_contribution(
    meta_founder, geno_probs, alternative_allele_prob, ind_mf
):
    """Split allele contributions by paternal and maternal metafounder source."""

    paternal_meta_founder = meta_founder[0]
    maternal_meta_founder = meta_founder[1]
    genotype1 = geno_probs[1, :]
    genotype2 = geno_probs[2, :]
    genotype3 = geno_probs[3, :]

    alternative_allele_prob[paternal_meta_founder] += 0.5 * genotype3 + genotype2
    alternative_allele_prob[maternal_meta_founder] += 0.5 * genotype3 + genotype1
    ind_mf[paternal_meta_founder] += 1
    ind_mf[maternal_meta_founder] += 1


def add_single_metafounder_contribution(
    meta_founder, geno_probs, alternative_allele_prob, ind_mf
):
    """Accumulate allele contributions for a founder with one metafounder source."""

    genotype1 = geno_probs[1, :]
    genotype2 = geno_probs[2, :]
    genotype3 = geno_probs[3, :]

    alternative_allele_prob[meta_founder] += 0.5 * (genotype2 + genotype1) + genotype3
    ind_mf[meta_founder] += 1


def update_maf_after_peeling(pedigree, peeling_info):
    """Updates the alternative allele frequency for each unknown parent group
    based on the mean genotype probabilities of the founders.
    Currently limits all alternative_allele_prob to be between 0.001 and 0.999.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function updates the pedigree.AAP attribute
    with the new alternative allele frequencies.
    """
    meta_founders = list(pedigree.AAP.keys())
    ind_mf = {k: 0 for k in meta_founders}
    alternative_allele_prob = {
        k: np.zeros(peeling_info.n_loci, dtype=np.float32) for k in meta_founders
    }
    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            add_founder_to_alt_allele_update(
                ind,
                peeling_info.get_geno_probs(ind.idn),
                alternative_allele_prob,
                ind_mf,
            )

    for mfx in meta_founders:
        current_aap = alternative_allele_prob[mfx]
        current_aap /= ind_mf[mfx]
        current_aap = np.maximum(np.minimum(current_aap, 0.999), 0.001)
        pedigree.AAP[mfx] = current_aap.astype(np.float32)

    maf_geno_cache = {}
    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            maf_geno = get_maf_genotypes_for_meta_founder(
                ind.MetaFounder, pedigree, peeling_info.n_loci, maf_geno_cache
            )
            peeling_info.anterior[ind.idn, :, :] = maf_geno


# Genotype and sequencing error-rate updates.


def update_penetrance(pedigree, peeling_info, args):
    """Updates the penetrance matrix for each individual.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :return: None. The function updates the peeling_info.penetrance attribute
    with the new genotype probabilities.
    """
    if args.est_geno_error_prob:
        peeling_info.geno_error = update_geno_error(pedigree, peeling_info)
    if args.est_seq_error_prob:
        peeling_info.seq_error = update_seq_error(pedigree, peeling_info)

    if peeling_info.is_x_chr:
        warnings.warn(
            "Updating error rates and minor allele frequencies for \
                X chromosomes is not well tested. \
                Recommend running without that option."
        )
    phase_founder = not args.no_phase_founder
    for ind in pedigree:
        x_chr_male_flag = peeling_info.is_x_chr and ind.sex == 0
        ind_penetrance = get_base_individual_penetrance(ind, peeling_info)

        if ind.phenotype is not None:
            ind_penetrance = updateGenoProbsFromPhenotype(
                ind_penetrance,
                ind.phenotype,
                pedigree.phenoPenetrance,
            )

        if ind.isGenotypedFounder() and phase_founder and ind.genotypes is not None:
            loci = get_het_midpoint(ind.genotypes)
            if loci is not None:
                error = peeling_info.geno_error[loci]
                if not x_chr_male_flag:
                    ind_penetrance[:, loci] = np.array(
                        [error / 3, error / 3, 1 - error, error / 3], dtype=np.float32
                    )

        peeling_info.penetrance[ind.idn, :, :] = ind_penetrance


def update_geno_error(pedigree, peeling_info):
    """Updates the genotype error rate for each locus using simple EM.
    Adds the expected number of errors that an individual has,
    marginalising over their current estimate of their genotype probabilities.
    We use a max value of 5% and a min value of .0001 percent to
    make sure the values are reasonable.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: the updated genotype error rates for each locus
    :rtype: 1D numpy array with length equal to the number of loci
    """
    counts = np.full(pedigree.nLoci, 1, dtype=np.float32)
    errors = np.full(pedigree.nLoci, 0.0001, dtype=np.float32)

    for ind in pedigree:
        update_geno_error_ind(
            counts, errors, ind.genotypes, peeling_info.get_geno_probs(ind.idn)
        )

    new_error = errors / counts
    new_error = np.maximum(np.minimum(new_error, 0.05), 0.0001)
    return new_error


@jit(nopython=True)
def update_geno_error_ind(counts, errors, genotypes, geno_probs):
    """Updates the genotype error rate at each locus with non-missing genotype.

    :param counts: vector of counts for each locus, initialized to 1
    :type counts: 1D numpy array with length equal to the number of loci
    :param errors: vector of errors for each locus, initialized to 0.001
    :type errors: 1D numpy array with length equal to the number of loci
    :param genotypes: observed genotypes for an individual collected
        via the geno_file input option.
    :type genotypes: 1D numpy array of Int8 with length n_loci
    :param geno_probs: genotype probabilities for each genotype state
        at each locus for the individual.
    :type geno_probs: 2D numpy array with shape 4 x n_loci
    :return: None. The function updates the counts and errors arrays in place.
    """
    for i, _ in enumerate(counts):
        if genotypes[i] != 9:  # Only include non-missing genotypes.
            counts[i] += 1
            genotype_prob0 = geno_probs[0, i]
            genotype_prob1 = geno_probs[1, i]
            genotype_prob2 = geno_probs[2, i]
            genotype_prob3 = geno_probs[3, i]
            if genotypes[i] == 0:
                errors[i] += genotype_prob1 + genotype_prob2 + genotype_prob3
            if genotypes[i] == 1:
                errors[i] += genotype_prob0 + genotype_prob3
            if genotypes[i] == 2:
                errors[i] += genotype_prob0 + genotype_prob1 + genotype_prob2


def update_seq_error(pedigree, peeling_info):
    """Updates the sequencing error rate at each locus homozygous states using simple EM.
    This update adds the expected number of errors that an individual
    has marginalizing over their current genotype probabilities.
    This only uses homozygous states; heterozygous states are
    ignored in both counts and errors.
    We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: the updated sequencing error rates for each locus
    :rtype: 1D numpy array with length equal to the number of loci
    """
    counts = np.full(pedigree.nLoci, 1, dtype=np.float32)
    errors = np.full(pedigree.nLoci, 0.001, dtype=np.float32)

    for ind in pedigree:
        if ind.reads is not None:
            update_seq_error_ind(
                counts,
                errors,
                ind.reads[0],
                ind.reads[1],
                peeling_info.get_geno_probs(ind.idn),
            )

    new_error = errors / counts
    new_error = np.maximum(np.minimum(new_error, 0.01), 0.0001)
    return new_error


@jit(nopython=True)
def update_seq_error_ind(counts, errors, ref_reads, alt_reads, geno_probs):
    """Updates the sequencing error rate for homozygous states at each locus with read data.

    :param counts: vector of counts for each locus, initialized to 1
    :type counts: 1D numpy array with length equal to the number of loci
    :param errors: vector of errors for each locus, initialized to 0.001
    :type errors: 1D numpy array with length equal to the number of loci
    :param ref_reads: the number of sequencing reads supporting the reference allele at each locus
    :type ref_reads: 1D numpy array of uint16 with length n_loci
    :param alt_reads: the number of sequencing reads supporting the alternative allele at each locus
    :type alt_reads: 1D numpy array of uint16 with length n_loci
    :param geno_probs: genotype probabilities for each genotype state
        at each locus for the individual.
    :type geno_probs: 2D numpy array with shape 4 x n_loci
    :return: None. The function updates the counts and errors arrays in place.
    """
    # Expected read errors are weighted by the probability of each homozygous genotype.
    for i, _ in enumerate(counts):
        genotype_prob0 = geno_probs[0, i]
        genotype_prob3 = geno_probs[3, i]
        counts[i] += (genotype_prob0 + genotype_prob3) * (alt_reads[i] + ref_reads[i])
        errors[i] += genotype_prob0 * alt_reads[i]
        errors[i] += genotype_prob3 * ref_reads[i]


def update_pheno_penetrance(pedigree, peeling_info):
    """Updates the phenotype penetrance matrix for each individual based on their phenotype.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. The function updates the pedigree.phenoPenetrance attribute
        with the new phenotype penetrance matrix.
    """
    # Based on Kinghorn (2003), "A Simple Method to Detect a Single Gene
    # that Determines a Categorical Trait with Incomplete Penetrance".
    rg_pheno = pedigree.phenoPenetrance.shape[1]
    denominator = np.full((4, pedigree.nLoci), 0, dtype=np.float32)
    contributions = np.full((4, rg_pheno), 0, dtype=np.float32)

    for ind in pedigree:
        if ind.phenotype is not None:
            update_pheno_penetrance_ind(
                denominator,
                contributions,
                rg_pheno,
                ind.phenotype,
                peeling_info.get_geno_probs(ind.idn),
            )

    for pheno in range(rg_pheno):
        pedigree.phenoPenetrance[:, pheno] = contributions[:, pheno] / denominator[:, 0]

    pedigree.phenoPenetrance = pedigree.phenoPenetrance / np.sum(
        pedigree.phenoPenetrance, 1, keepdims=True
    )


def update_pheno_penetrance_ind(
    denominator, contributions, rg_pheno, phenotype, geno_probs
):
    """Updates the phenotype penetrance matrix for an individual based on
    their phenotype and genotype probabilities.

    :param denominator: Sums the genotype probabilities across individuals
        with phenotype data, initialised to 0.
    :type denominator: 2D numpy array with shape 4 x n_loci
    :param contributions: matrix of contributions for each phenotype
        and genotype state, initialized to 0
    :type contributions: 2D numpy array with shape nPhenotype categories x 4
    :param rg_pheno: number of phenotype categories
    :type rg_pheno: int
    :param phenotype: the phenotype of the individual
    :type phenotype: int
    :param geno_probs: genotype probabilities for each genotype state
        at each locus for the individual.
    :type geno_probs: 2D numpy array with shape 4 x n_loci
    :return: None. The function updates the counts and contributions arrays in place.
    """
    # Multiple phenotype records for an individual contribute as repeated counts.
    geno_probs_first_locus = geno_probs[:, 0]
    for pheno in phenotype:
        pheno = int(pheno)
        if 0 <= pheno < rg_pheno:
            denominator += geno_probs
            contributions[:, pheno] += geno_probs_first_locus
