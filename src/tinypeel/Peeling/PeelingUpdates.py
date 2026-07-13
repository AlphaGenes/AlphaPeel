from numba import jit
import numpy as np

from ..tinyhouse import ProbMath

from . import PeelingInfo

import warnings

#########################################################################################
# In this module we will update 3 things:                                               #
# 1) Our estimate for the MAF (both prior to peeling and after each peeling cycle) .    #
# 2) Our estimate of the locus specific (sequencing) error rate.                        #
# 3) Our estimate of the locus specific recombination rate.                             #
#########################################################################################


# Estimating the alternative allele frequency. This update is done by using an iterative approach
# which maximizes the likelihood of the observed genotypes conditional on them having been
# generated from hardy-weinberg equilibrium with a fixed maf value. To speed up, we use
# Newton style updates to estimate the alternative allele frequency.
# There is math on how to do this... somewhere?


def updateMaf(pedigree, peelingInfo):
    """Estimates the alternative allele frequency at all loci (i.e markers).

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. The function updates the pedigree.AAP attribute with the new alternative allele frequencies.
    """
    if peelingInfo.isXChr:
        warnings.warn(
            "Updating error rates and alternative allele frequencies for X chromosomes are not well test and will break in interesting ways. Recommend running without that option."
        )
    MF = list(pedigree.AAP.keys())
    for mfx in MF:
        AAP = pedigree.AAP[mfx]
        for i in range(peelingInfo.nLoci):
            AAP[i] = newtonMafUpdates(peelingInfo, AAP, i)
        pedigree.AAP[mfx] = AAP.astype(np.float32)

    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            AAP = {
                k: np.zeros(peelingInfo.nLoci, dtype=np.float32)
                for k in ind.MetaFounder
            }
            for mfx in ind.MetaFounder:
                AAP[mfx] = pedigree.AAP[mfx]
            if len(ind.MetaFounder) == 2:
                mafGeno = ProbMath.getGenotypesFromMultiMaf(AAP)
            else:
                mafGeno = ProbMath.getGenotypesFromMaf(AAP[mfx])
            peelingInfo.anterior[ind.idn, :, :] = mafGeno


def newtonMafUpdates(peelingInfo, AAP, index):
    """Iterative approximation for the prior alternative allele frequency.
    Currently limits all AAP to be between 0.001 and 0.999.

    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param AAP: starting alternative allele frequencies, default 0.5
    :type AAP: 1D numpy array with length equal to the number of loci
    :param index: the marker index for which to update the alternative allele frequency
    :type index: int
    :return: the updated alternative allele frequency for the given marker index
    :rtype: float
    """

    if AAP[index] < 0.001:
        maf = 0.001
    elif AAP[index] > 0.999:
        maf = 0.999
    else:
        maf = AAP[index]

    iters = 5
    converged = False
    while not converged:
        maf_old = maf
        delta = getNewtonUpdate(maf_old, peelingInfo, index)
        maf = maf_old + delta
        if maf < 0.001:
            maf = 0.001
        if maf > 0.999:
            maf = 0.999
        if abs(maf - maf_old) < 0.0001:
            converged = True
        iters -= 1
        if iters < 0:
            converged = True
    return maf


@jit(nopython=True)
def getNewtonUpdate(p, peelingInfo, index):
    """Calculates the alternative allele frequency using Newton's method of optimisation.

    :param p: the current alternative allele frequency estimate
    :type p: float
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param index: the marker index for which to update the alternative allele frequency
    :type index: int
    :return: ratio of the first and second derivatives of the log likelihood function to be added to the current alternative allele frequency estimate.
    :rtype: float
    """
    # Log liklihood's first + second derivitives
    LLp = 0
    LLpp = 0

    # I want to add priors. Should be 1 individual of each of the four states.
    LLp, LLpp = addIndividualToUpdate(
        np.array([1, 0, 0, 0], dtype=np.float32), p, LLp, LLpp
    )
    LLp, LLpp = addIndividualToUpdate(
        np.array([0, 1, 0, 0], dtype=np.float32), p, LLp, LLpp
    )
    LLp, LLpp = addIndividualToUpdate(
        np.array([0, 0, 1, 0], dtype=np.float32), p, LLp, LLpp
    )
    LLp, LLpp = addIndividualToUpdate(
        np.array([0, 0, 0, 1], dtype=np.float32), p, LLp, LLpp
    )
    for i in range(peelingInfo.nInd):
        if peelingInfo.genotyped[i, index]:
            d = peelingInfo.penetrance[i, :, index]
            LLp, LLpp = addIndividualToUpdate(d, p, LLp, LLpp)
    if LLp == 0 or LLpp == 0:
        return 0  # Could be a case where no one has data.
    return -LLp / LLpp


@jit(nopython=True)
def addIndividualToUpdate(d, p, LLp, LLpp):
    """Adds each available genotype data to first and second derivatives of the log likelihood.

    :param d: the penetrance term for genotyped individuals as genotype probabilities
    :type d: 1D numpy array with length 4 (i.e [p(AA), p(aA), p(Aa), p(aa)])
    :param p: the current alternative allele frequency estimate
    :type p: float
    :param LLp: the first derivative of the log likelihood
    :type LLp: float
    :param LLpp: the second derivative of the log likelihood
    :type LLpp: float
    :return: updated first and second derivatives of the log likelihood
    :rtype: tuple(float, float)
    """
    d0 = d[0]
    d1 = d[1] + d[2]
    d2 = d[3]

    f = d0 * (1 - p) ** 2 + d1 * p * (1 - p) + d2 * p**2
    fp = (d1 - 2 * d0) + 2 * p * (d0 + d2 - d1)
    fpp = 2 * (d0 + d2 - d1)

    LLp += fp / f
    LLpp += fpp / f - (fp / f) ** 2

    return LLp, LLpp


def updateMafAfterPeeling(pedigree, peelingInfo):
    """Updates the alternative allele frequency for each unknown parent group based on the mean genotype probabilities of the founders.
    Currently limits all AAP to be between 0.001 and 0.999.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. The function updates the pedigree.AAP attribute with the new alternative allele frequencies.
    """
    MF = list(pedigree.AAP.keys())
    indMF = {k: 0 for k in MF}
    AAP = {k: np.zeros(peelingInfo.nLoci, dtype=np.float32) for k in MF}
    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            ind_genotype = peelingInfo.getGenoProbs(ind.idn)
            if len(ind.MetaFounder) == 2:
                AAP[ind.MetaFounder[0]] += 0.5 * ind_genotype[3, :] + ind_genotype[2, :]
                AAP[ind.MetaFounder[1]] += 0.5 * ind_genotype[3, :] + ind_genotype[1, :]
                indMF[ind.MetaFounder[0]] += 1
                indMF[ind.MetaFounder[1]] += 1
            else:
                AAP[ind.MetaFounder[0]] += (
                    0.5 * (ind_genotype[2, :] + ind_genotype[1, :]) + ind_genotype[3, :]
                )
                indMF[ind.MetaFounder[0]] += 1

    for mfx in MF:
        for i in range(peelingInfo.nLoci):
            AAP[mfx][i] = AAP[mfx][i] / indMF[mfx]
            if AAP[mfx][i] < 0.001:
                AAP[mfx][i] = 0.001
            elif AAP[mfx][i] > 0.999:
                AAP[mfx][i] = 0.999
        pedigree.AAP[mfx] = AAP[mfx].astype(np.float32)

    for ind in pedigree:
        if ind.MetaFounder is not None and ind.isFounder():
            AAP = {
                k: np.zeros(peelingInfo.nLoci, dtype=np.float32)
                for k in ind.MetaFounder
            }
            for mfx in ind.MetaFounder:
                AAP[mfx] = pedigree.AAP[mfx]
            if len(ind.MetaFounder) == 2:
                mafGeno = ProbMath.getGenotypesFromMultiMaf(AAP)
            else:
                mafGeno = ProbMath.getGenotypesFromMaf(AAP[mfx])
            peelingInfo.anterior[ind.idn, :, :] = mafGeno


#
# NOTE: The following code updates the genotype and sequencing error rates.
#


def updatePenetrance(pedigree, peelingInfo, args):
    """Updates the penetrance matrix for each individual.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :return: None. The function updates the peelingInfo.penetrance attribute with the new genotype probabilities.
    """
    if args.est_geno_error_prob:
        peelingInfo.genoError = updateGenoError(pedigree, peelingInfo)
    if args.est_seq_error_prob:
        peelingInfo.seqError = updateSeqError(pedigree, peelingInfo)

    if peelingInfo.isXChr:
        warnings.warn(
            "Updating error rates and minor allele frequencies for X chromosomes are not well test and will break in interesting ways. Recommend running without that option."
        )
    phaseFounder = not args.no_phase_founder
    for ind in pedigree:
        XChrMaleFlag = (
            peelingInfo.isXChr and ind.sex == 0
        )  # This is the X chromosome and the individual is male.
        peelingInfo.penetrance[ind.idn, :, :] = ProbMath.getGenotypeProbabilities(
            peelingInfo.nLoci,
            ind.genotypes,
            ind.reads,
            peelingInfo.genoError,
            peelingInfo.seqError,
            XChrMaleFlag,
        )

        if ind.phenotype is not None:
            peelingInfo.penetrance[
                ind.idn, :, :
            ] = ProbMath.updateGenoProbsFromPhenotype(
                peelingInfo.penetrance[ind.idn, :, :],
                ind.phenotype,
                pedigree.phenoPenetrance,
            )

        if ind.isGenotypedFounder() and phaseFounder and ind.genotypes is not None:
            loci = PeelingInfo.getHetMidpoint(ind.genotypes)
            if loci is not None:
                error = peelingInfo.genoError[loci]
                if not XChrMaleFlag:
                    peelingInfo.penetrance[ind.idn, :, loci] = np.array(
                        [error / 3, error / 3, 1 - error, error / 3], dtype=np.float32
                    )


def updateGenoError(pedigree, peelingInfo):
    """Updates the genotype error rate for each locus using simple EM.
    Adds the expected number of errors that an individual has, marginalising over their current estimate of their genotype probabilities.
    We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: the updated genotype error rates for each locus
    :rtype: 1D numpy array with length equal to the number of loci
    """
    counts = np.full(pedigree.nLoci, 1, dtype=np.float32)
    errors = np.full(pedigree.nLoci, 0.0001, dtype=np.float32)

    for ind in pedigree:
        updateGenoError_ind(
            counts, errors, ind.genotypes, peelingInfo.getGenoProbs(ind.idn)
        )

    newError = errors / counts
    newError = np.maximum(np.minimum(newError, 0.05), 0.0001)
    return newError


@jit(nopython=True)
def updateGenoError_ind(counts, errors, genotypes, genoProbs):
    """Updates the genotype error rate at each locus with non-missing genotype.

    :param counts: vector of counts for each locus, initialized to 1
    :type counts: 1D numpy array with length equal to the number of loci
    :param errors: vector of errors for each locus, initialized to 0.001
    :type errors: 1D numpy array with length equal to the number of loci
    :param genotypes: observed genotypes for an individual collected via the geno_file input option.
    :type genotypes: 1D numpy array of Int8 with length nLoci
    :param genoProbs: genotype probabilities for each genotype state at each locus for the individual.
    :type genoProbs: 2D numpy array with shape 4 x nLoci
    :return: None. The function updates the counts and errors arrays in place.
    """
    for i in range(len(counts)):
        if genotypes[i] != 9:  # Only include non-missing genotypes.
            counts[i] += 1
            if genotypes[i] == 0:
                errors[i] += genoProbs[1, i] + genoProbs[2, i] + genoProbs[3, i]
            if genotypes[i] == 1:
                errors[i] += genoProbs[0, i] + genoProbs[3, i]
            if genotypes[i] == 2:
                errors[i] += genoProbs[0, i] + genoProbs[1, i] + genoProbs[2, i]


def updateSeqError(pedigree, peelingInfo):
    """Updates the sequencing error rate at each locus homozygous states using simple EM.
    This update adds the expected number of errors that an individual has marginalizing over their current genotype probabilities.
    This only uses the homozygotic states, heterozygotic states are ignored (in both the counts + errors terms).
    We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: the updated sequencing error rates for each locus
    :rtype: 1D numpy array with length equal to the number of loci
    """
    counts = np.full(pedigree.nLoci, 1, dtype=np.float32)
    errors = np.full(pedigree.nLoci, 0.001, dtype=np.float32)

    for ind in pedigree:
        if ind.reads is not None:
            updateSeqError_ind(
                counts,
                errors,
                ind.reads[0],
                ind.reads[1],
                peelingInfo.getGenoProbs(ind.idn),
            )

    newError = errors / counts
    newError = np.maximum(np.minimum(newError, 0.01), 0.0001)
    return newError


@jit(nopython=True)
def updateSeqError_ind(counts, errors, refReads, altReads, genoProbs):
    """Updates the sequencing error rate for homozygotic states at each locus with non-missing genotype.

    :param counts: vector of counts for each locus, initialized to 1
    :type counts: 1D numpy array with length equal to the number of loci
    :param errors: vector of errors for each locus, initialized to 0.001
    :type errors: 1D numpy array with length equal to the number of loci
    :param refReads: the number of sequencing reads supporting the reference allele at each locus
    :type refReads: 1D numpy array of int64 with length nLoci
    :param altReads: the number of sequencing reads supporting the alternative allele at each locus
    :type altReads: 1D numpy array of int64 with length nLoci
    :param genoProbs: genotype probabilities for each genotype state at each locus for the individual.
    :type genoProbs: 2D numpy array with shape 4 x nLoci
    :return: None. The function updates the counts and errors arrays in place.
    """
    # Errors occur when genotype is 0 and an alternative allele happens.
    # Errors occur when genotype is 2 (coded as 3) and a reference allele happens.
    # Number of observations is number of reads * probability the individual is homozygous.
    for i in range(len(counts)):
        counts[i] += (genoProbs[0, i] + genoProbs[3, i]) * (altReads[i] + refReads[i])
        errors[i] += genoProbs[0, i] * altReads[i]
        errors[i] += genoProbs[3, i] * refReads[i]


def updatePhenoPenetrance(pedigree, peelingInfo):
    """Updates the phenotype penetrance matrix for each individual based on their phenotype.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. The function updates the pedigree.phenoPenetrance attribute with the new phenotype penetrance matrix.
    """
    # Credit to Kinghorn (2003) A SIMPLE METHOD TO DETECT A SINGLE GENE THAT DETERMINES ACATEGORICAL TRAIT WITH INCOMPLETE PENETRANCE
    rgPheno = pedigree.phenoPenetrance.shape[1]  # Range of phenotype values
    denominator = np.full(
        (4, pedigree.nLoci), 0, dtype=np.float32
    )  # Sum of the genotypes across individuals with any phenotype data
    contributions = np.full((4, rgPheno), 0, dtype=np.float32)

    for ind in pedigree:
        if ind.phenotype is not None:
            updatePhenoPenetrance_ind(
                denominator,
                contributions,
                rgPheno,
                ind.phenotype,
                peelingInfo.getGenoProbs(ind.idn),
            )

    for pheno in range(rgPheno):
        pedigree.phenoPenetrance[:, pheno] = contributions[:, pheno] / denominator[:, 0]

    # Normalize the contributions to get the penetrance matrix.
    pedigree.phenoPenetrance = pedigree.phenoPenetrance / np.sum(
        pedigree.phenoPenetrance, 1, keepdims=True
    )


def updatePhenoPenetrance_ind(
    denominator, contributions, rgPheno, phenotype, genoProbs
):
    """Updates the phenotype penetrance matrix for an individual based on their phenotype and genotype probabilities.

    :param denominator: Sums the genotype probabilities across individuals with phenotype data, initialised to 0.
    :type denominator: 2D numpy array with shape 4 x nLoci
    :param contributions: matrix of contributions for each phenotype and genotype state, initialized to 0
    :type contributions: 2D numpy array with shape nPhenotype categories x 4
    :param rgPheno: number of phenotype categories
    :type rgPheno: int
    :param phenotype: the phenotype of the individual
    :type phenotype: int
    :param genoProbs: genotype probabilities for each genotype state at each locus for the individual.
    :type genoProbs: 2D numpy array with shape 4 x nLoci
    :return: None. The function updates the counts and contributions arrays in place.
    """
    # For now, assuming only single locus genotype input
    # Handles multiple phenotype record as another count

    for pheno in phenotype:
        pheno = int(pheno)
        if 0 <= pheno < rgPheno:
            denominator += genoProbs
            contributions[:, pheno] += genoProbs[:, 0]
