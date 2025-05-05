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
    """
    This function estimates the alternative allele frequency at all loci (i.e markers).
    Updates the alternative allele frequency for each unknown parent group, then the anterior term for each founder.
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
    """
    This function gives an iterative approximation for the alternative allele frequency.
    First the starting alternative allele frequency is restricted to be between 0.01 and 0.99.
    Then the iteration is set up for Newton's method. This finishes at convergence or 5 iterations.
    Each iteration collects delta via the getNewtonUpdate function, adds delta to the current alternative alllele frequency,
    and checks for convergence.
    The function returns the final alternative allele frequency.
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
    """
    This function calculates the Newton update for the alternative allele frequency.
    Considers uncertainty in the available genotype data through the penetrance.
    Returns delta
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
    """
    This function adds each available genotype data to first and second derivatives of the log likelihood.
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
    """
    This function updates the alternative allele frequency for each unknown parent group based on the mean genotype probabilities of the founders.
    If triggered, this function will be called after each peeling cycle.
    Currently restricts alternative allele frequencies from 0.01 to 0.99.
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
    """
    This function updates the penetrance matrix for each individual.
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
    """
    This function updates the genotype error rate for each locus using simple EM.
    Adds the expected number of errors that an individual has, marginalising over their current estimate of their genotype probabilities.
    We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable.
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
    """
    This function updates the genotype error rate at each locus with non-missing genotype.
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
    """
    The following is a simple EM update for the genotyping error rate at each locus.
    This update adds the expected number of errors that an individual has marginalizing over their current genotype probabilities.
    This only uses the homozygotic states, heterozygotic states are ignored (in both the counts + errors terms).
    We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable
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
    """
    This function updates the sequencing error rate at each locus with non-missing genotype.
    This is completed only for the homozygotic states.
    """
    # Errors occur when genotype is 0 and an alternative allele happens.
    # Errors occur when genotype is 2 (coded as 3) and a reference allele happens.
    # Number of observations is number of reads * probability the individual is homozygous.
    for i in range(len(counts)):
        counts[i] += (genoProbs[0, i] + genoProbs[3, i]) * (altReads[i] + refReads[i])
        errors[i] += genoProbs[0, i] * altReads[i]
        errors[i] += genoProbs[3, i] * refReads[i]


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
    """
    This function estimates the distance of a locus from a break in the haplotype due to recombination.
    """
    breakPoints = np.floor(np.linspace(loc1, loc2, nBreaks)).astype(dtype=np.int64)
    distances = [
        getDistance(breakPoints[i], breakPoints[i + 1], peelingInfo)
        for i in range(nBreaks - 1)
    ]
    print(nBreaks, ":", distances)
    return sum(distances)


def getDistance(loc1, loc2, peelingInfo):
    """
    This function estimates the distance between two loci based on the recombination rate.
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
    """
    This function calculates the total probability of the grand maternal allele being transmitted from the father in patSeg.
    """
    seg = peelingInfo.segregation[:, :, loc]
    sumSeg = np.sum(seg, 1)
    patSeg = (seg[:, 2] + seg[:, 3]) / sumSeg
    # matSeg = (seg[:,1] + seg[:,3])/sumSeg
    return patSeg


def haldane(difference):
    """
    This function converts the recombination rate to a distance using Haldane's formula.
    Haldane's formula considers the fixation probability of a beneficial allele in a population.
    """
    return -np.log(1.0 - 2.0 * difference) / 2.0
