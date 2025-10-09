from numba import jit
import numpy as np

from ..tinyhouse import InputOutput
from ..tinyhouse import ProbMath

from . import PeelingInfo

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
    if peelingInfo.isSexChrom:
        print(
            "Updating error rates and alternative allele frequencies for sex chromosomes are not well test and will break in interesting ways. Recommend running without that option."
        )
    MF = list(pedigree.AAP.keys())
    for mfx in MF:
        AAP = pedigree.AAP[mfx]
        for i in range(peelingInfo.nLoci):
            AAP[i] = newtonMafUpdates(peelingInfo, AAP, i)

        mafGeno = ProbMath.getGenotypesFromMaf(AAP)
        for ind in pedigree:
            if ind.MetaFounder == mfx and ind.isFounder():
                peelingInfo.anterior[ind.idn, :, :] = mafGeno
        pedigree.AAP[mfx] = AAP.astype(np.float32)


def newtonMafUpdates(peelingInfo, AAP, index):
    """Iterative approximation for the prior alternative allele frequency.
    Currently limits all AAP to be between 0.01 and 0.99.

    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :param AAP: starting alternative allele frequencies, default 0.5
    :type AAP: 1D numpy array with length equal to the number of loci
    :param index: the marker index for which to update the alternative allele frequency
    :type index: int
    :return: the updated alternative allele frequency for the given marker index
    :rtype: float
    """

    if AAP[index] < 0.01:
        maf = 0.01
    elif AAP[index] > 0.99:
        maf = 0.99
    else:
        maf = AAP[index]

    iters = 5
    converged = False
    while not converged:
        maf_old = maf
        delta = getNewtonUpdate(maf_old, peelingInfo, index)
        maf = maf_old + delta
        if maf < 0.01:
            maf = 0.01
        if maf > 0.99:
            maf = 0.99
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
    :type d: 1D numpy array with length 4 (i.e [p(AA), p(Aa), p(aA), p(aa)])
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
    Currently limits all AAP to be between 0.01 and 0.99.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: None. The function updates the pedigree.AAP attribute with the new alternative allele frequencies.
    """
    MF = list(pedigree.AAP.keys())
    for mfx in MF:
        AAP = pedigree.AAP[mfx]
        sumOfGenotypes = np.full((4, peelingInfo.nLoci), 0, dtype=np.float32)
        n = 0
        for ind in pedigree:
            if ind.MetaFounder == mfx and ind.isFounder():
                ind_genotype = peelingInfo.getGenoProbs(ind.idn)
                sumOfGenotypes += ind_genotype
                n += 1
        for i in range(peelingInfo.nLoci):
            marker_Aa = sumOfGenotypes[1, i]
            marker_aA = sumOfGenotypes[2, i]
            marker_AA = sumOfGenotypes[3, i]
            AAP[i] = (0.5 * (marker_Aa + marker_aA) + marker_AA) / n
            if AAP[i] < 0.01:
                AAP[i] = 0.01
            elif AAP[i] > 0.99:
                AAP[i] = 0.99

        mafGeno = ProbMath.getGenotypesFromMaf(AAP)
        for ind in pedigree:
            if ind.MetaFounder == mfx and ind.isFounder():
                peelingInfo.anterior[ind.idn, :, :] = mafGeno
        pedigree.AAP[mfx] = AAP.astype(np.float32)


# Commenting out the following code. This was used to do updates via grid search.
# @jit(nopython = True)
# def mafLoglikelihood(peelingInfo, maf, index):
#     score = 0
#     maf_squared = maf**2
#     maf_one_minus_maf = maf*(1-maf)
#     one_minus_maf_squared = (1-maf)**2

#     for i in range(peelingInfo.nInd):
#         if peelingInfo.genotyped[i, index] :
#         # if True :
#             genoProbs = peelingInfo.penetrance[i,:,index]

#             prob = 0
#             prob += one_minus_maf_squared*genoProbs[0]
#             prob += maf_one_minus_maf*genoProbs[1]
#             prob += maf_one_minus_maf*genoProbs[2]
#             prob += maf_squared*genoProbs[3]
#             score += math.log(prob)
#     return score

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

    if peelingInfo.isSexChrom:
        print(
            "Updating error rates and minor allele frequencies for sex chromosomes are not well test and will break in interesting ways. Recommend running without that option."
        )
    for ind in pedigree:
        sexChromFlag = (
            peelingInfo.isSexChrom and ind.sex == 0
        )  # This is the sex chromosome and the individual is male.
        peelingInfo.penetrance[ind.idn, :, :] = ProbMath.getGenotypeProbabilities(
            peelingInfo.nLoci,
            ind.genotypes,
            ind.reads,
            peelingInfo.genoError,
            peelingInfo.seqError,
            sexChromFlag,
        )

        if ind.phenotype is not None:
            peelingInfo.penetrance[
                ind.idn, :, :
            ] = ProbMath.updateGenoProbsFromPhenotype(
                peelingInfo.penetrance[ind.idn, :, :],
                ind.phenotype,
                pedigree.phenoPenetrance,
            )

        if (
            ind.isGenotypedFounder()
            and (not InputOutput.args.no_phase_founder)
            and ind.genotypes is not None
        ):
            loci = PeelingInfo.getHetMidpoint(ind.genotypes)
            if loci is not None:
                e = peelingInfo.genoError[loci]
                peelingInfo.penetrance[ind.idn, :, loci] = np.array(
                    [e / 3, e / 3, 1 - e, e / 3], dtype=np.float32
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
    rgPheno = len(pedigree.phenoPenetrance[0, :])  # Range of phenotype values
    denominator =  np.full((4, pedigree.nLoci), 0, dtype = np.float32) # Sum of the genotypes across individuals with any phenotype data
    #counts = np.full(rgPheno, 0, dtype=np.float32)
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

    #mask = counts > 0
    #pedigree.phenoPenetrance[:, mask] = contributions[:, mask] / counts[mask]
    
    for pheno in range(rgPheno):
        pedigree.phenoPenetrance[:, pheno] = contributions[:, pheno] / denominator[:,0]

    
    # Normalize the contributions to get the penetrance matrix.
    pedigree.phenoPenetrance = pedigree.phenoPenetrance / np.sum(
        pedigree.phenoPenetrance, 1, keepdims=True
    )
    

def updatePhenoPenetrance_ind(denominator, contributions, rgPheno, phenotype, genoProbs):
    """Updates the phenotype penetrance matrix for an individual based on their phenotype and genotype probabilities.

    :param denominator: 
    :type denominator: 
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


#
# NOTE: The following code is designed to estimate the recombination rate between markers.
# This is currently broken and does not give good estimates. Strongly recommend not using it, and the associated options have been disabled in tinyPeel.py
#


# def updateSeg(peelingInfo):
#     # Idea: Split the chromosome into ~10 blocks with slightly different starts/stops. Estimate at each one and then sum.
#     # Not sure if this is the best way, but probably will work okay.

#     distances = [
#         estDistanceFromBreaks(0, peelingInfo.nLoci - 1, nBreaks, peelingInfo)
#         for nBreaks in [5, 10, 20]
#     ]
#     distance = np.mean(distances[0:2])
#     print(distances, ":", distance)

#     setupTransmission(distance, peelingInfo)


def estDistanceFromBreaks(loc1, loc2, nBreaks, peelingInfo):
    """Estimates the distance of a locus from a break in the haplotype due to recombination.

    :param loc1: first locus number
    :type loc1: float
    :param loc2: second locus number
    :type loc2: float
    :param nBreaks: number of breaks to estimate the distance between the two loci.
    :type nBreaks: int
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: the estimated distance between the two loci based on the recombination rate.
    :rtype: float
    """
    breakPoints = np.floor(np.linspace(loc1, loc2, nBreaks)).astype(dtype=np.int64)
    distances = [
        getDistance(breakPoints[i], breakPoints[i + 1], peelingInfo)
        for i in range(nBreaks - 1)
    ]
    print(nBreaks, ":", distances)
    return sum(distances)


def getDistance(loc1, loc2, peelingInfo):
    """Estimates the distance between two loci based on the recombination rate.

    :param loc1: first locus number
    :type loc1: float
    :param loc2: second locus number
    :type loc2: float
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: the estimated recombination rate between two loci in distance units.
    :rtype: float
    """
    patSeg1 = getSumSeg(loc1, peelingInfo)
    patSeg2 = getSumSeg(loc2, peelingInfo)

    # patValid = ((patSeg1 > 0.99) | (patSeg1 < 0.01)) & (
    #     (patSeg2 > 0.99) | (patSeg2 < 0.01)
    # )

    # matValid = ((matSeg1 > .99) | (matSeg1 < .01)) & ((matSeg2 > .99) | (matSeg2 < .01))

    # patRecomb = np.mean(np.abs(patSeg1[patValid] - patSeg2[patValid]))
    # matRecomb = np.mean(np.abs(matSeg1[matValid] - matSeg2[matValid]))

    # recomb = (patRecomb + matRecomb)/2

    # difference = np.abs(np.round(patSeg1[patValid]) - np.round(patSeg2[patValid]))
    # return np.mean(difference)

    entropy1 = 1 - (-patSeg1 * np.log2(patSeg1) - (1 - patSeg1) * np.log2(1 - patSeg1))
    entropy2 = 1 - (-patSeg2 * np.log2(patSeg2) - (1 - patSeg2) * np.log2(1 - patSeg2))

    difference = patSeg1 * (1 - patSeg2) + (1 - patSeg1) * patSeg2
    est = np.sum(entropy1 * entropy2 * difference) / np.sum(entropy1 * entropy2)
    # return est
    return haldane(est)


def getSumSeg(loc, peelingInfo):
    """Calculates the total probability of the grand maternal allele being transmitted from the father in patSeg.

    :param loc: locus/marker number
    :type loc: int
    :param peelingInfo: Peeling information container.
    :type peelingInfo: class:`PeelingInfo.jit_peelingInformation`
    :return: the probability of the grand maternal allele being transmitted from the father at the given locus.
    :rtype: 1D numpy array with length equal to the number of individuals
    """
    seg = peelingInfo.segregation[:, :, loc]
    sumSeg = np.sum(seg, 1)
    patSeg = (seg[:, 2] + seg[:, 3]) / sumSeg
    # matSeg = (seg[:,1] + seg[:,3])/sumSeg
    return patSeg


def haldane(difference):
    """Converts the recombination rate to a distance using Haldane's formula.
    Haldane's formula considers the fixation probability of a beneficial allele in a population.

    :param difference: estimated recombination rate between two loci.
    :type difference: float
    :return: Haldane distance, which is the distance between two loci based on the recombination rate.
    :rtype: float
    """
    return -np.log(1.0 - 2.0 * difference) / 2.0
