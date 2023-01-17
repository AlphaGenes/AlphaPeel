import concurrent.futures
from numba import jit, float32, int8, int64, optional, boolean
from numba.experimental import jitclass
import numpy as np
from collections import OrderedDict

from ..tinyhouse import InputOutput
from ..tinyhouse import ProbMath
from ..tinyhouse import HaplotypeOperations

from . import PeelingInfo

import math

#########################################################################
### In this module we will update 3 things:                           ###  
### 1) Our estimate for the MAF.                                      ### 
### 2) Our estimate of the locus specific (sequencing) error rate.    ###
### 3) Our estimate of the locus specific recombination rate.         ###
#########################################################################


# Estimating the minor allele frequency. This update is done by using an iterative approach 
# which minimizes the likelihood of the observed genotypes conditional on them having been 
# generated from hardy-weinberg equilibrium with a fixed maf value. To speed up, we use
# Newton style updates to re-estimate the minor allele frequency.
# There is math on how to do this... somewhere?

def updateMaf(pedigree, peelingInfo):
    if peelingInfo.isSexChrom:
        print("Updating error rates and minor allele frequencies for sex chromosomes are not well test and will break in interesting ways. Recommend running without that option.")

    maf = np.full(peelingInfo.nLoci, .5, dtype = np.float32)
    for i in range(peelingInfo.nLoci):
        maf[i] = newtonMafUpdates(peelingInfo, i)

    mafGeno = ProbMath.getGenotypesFromMaf(maf)
    for ind in pedigree:
        if ind.isFounder():
            peelingInfo.anterior[ind.idn,:,:] = mafGeno
    peelingInfo.maf = maf.astype(np.float32)

def newtonMafUpdates(peelingInfo, index) :
    # This function gives an iterative approximation for the minor allele frequency. It uses a maximum of 5 iterations.
    maf = 0.5
    maf_old = 0.5
    iters = 5
    converged = False
    while not converged:
        delta = getNewtonUpdate(maf_old, peelingInfo, index)
        maf = maf_old + delta
        if maf < 0.01:
            maf = 0.01
        if maf > .99:
            maf = .99
        if abs(maf - maf_old) < 0.0001:
            converged = True
        iters -= 1
        if iters < 0:
            converged = True
    return maf

@jit(nopython=True)
def getNewtonUpdate(p, peelingInfo, index) :
    #Log likelihood of the liklihood's first + second derivitives
    LLp = 0
    LLpp = 0

    #I want to add priors. Should be 1 individual of each of the four states.
    LLp, LLpp = addIndividualToUpdate(np.array([1, 0, 0, 0], dtype = np.float32), p, LLp, LLpp)
    LLp, LLpp = addIndividualToUpdate(np.array([0, 1, 0, 0], dtype = np.float32), p, LLp, LLpp)
    LLp, LLpp = addIndividualToUpdate(np.array([0, 0, 1, 0], dtype = np.float32), p, LLp, LLpp)
    LLp, LLpp = addIndividualToUpdate(np.array([0, 0, 0, 1], dtype = np.float32), p, LLp, LLpp)
    for i in range(peelingInfo.nInd):
        if peelingInfo.genotyped[i, index] :
            d = peelingInfo.penetrance[i,:,index]
            LLp, LLpp = addIndividualToUpdate(d, p, LLp, LLpp)
    if LLp == 0 or LLpp == 0: return 0 #Could be a case where no one has data.
    return -LLp/LLpp

@jit(nopython = True)
def addIndividualToUpdate(d, p, LLp, LLpp):
    d0 = d[0]
    d1 = d[1] + d[2]
    d2 = d[3]

    f = d0*(1-p)**2 + d1*p*(1-p) + d2*p**2
    fp = (d1 - 2*d0) + 2*p*(d0 + d2 - d1)
    fpp = 2*(d0 + d2 - d1)
    
    LLp += fp/f
    LLpp += fpp/f - (fp/f)**2

    return LLp, LLpp

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

###
### NOTE: The following code updates the genotype and sequencing error rates.
###

def updatePenetrance(pedigree, peelingInfo):
    peelingInfo.genoError = updateGenoError(pedigree, peelingInfo)
    peelingInfo.seqError = updateSeqError(pedigree, peelingInfo)

    if peelingInfo.isSexChrom:
        print("Updating error rates and minor allele frequencies for sex chromosomes are not well test and will break in interesting ways. Recommend running without that option.")
    for ind in pedigree:
        sexChromFlag = peelingInfo.isSexChrom and ind.sex == 0 #This is the sex chromosome and the individual is male.
        peelingInfo.penetrance[ind.idn,:,:] = ProbMath.getGenotypeProbabilities(peelingInfo.nLoci, ind.genotypes, ind.reads, peelingInfo.genoError, peelingInfo.seqError, sexChromFlag)

        if ind.isGenotypedFounder() and (not InputOutput.args.nophasefounders) and ind.genotypes is not None:
            loci = PeelingInfo.getHetMidpoint(ind.genotypes)
            if loci is not None:
                e = peelingInfo.genoError[loci]
                peelingInfo.penetrance[ind.idn,:,loci] = np.array([e/3, e/3, 1-e, e/3], dtype = np.float32)

def updateGenoError(pedigree, peelingInfo) :
    # The following is a simple EM update for the genotyping error rate at each locus.
    # This update adds the expected number of errors that an individual marginalizing over their current estimate of their genotype probabilities.
    # We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable

    counts = np.full(pedigree.nLoci, 1, dtype = np.float32)
    errors = np.full(pedigree.nLoci, 0.01, dtype = np.float32)

    for ind in pedigree:
        updateGenoError_ind(counts, errors, ind.genotypes, peelingInfo.getGenoProbs(ind.idn))
    
    newError = errors/counts 
    newError = np.maximum(np.minimum(newError, .05), .0001)
    return newError

@jit(nopython=True)
def updateGenoError_ind(counts, errors, genotypes, genoProbs):
    for i in range(len(counts)) :
        if genotypes[i] != 9: # Only include non-missing genotypes.
            counts[i] += 1
            if genotypes[i] == 0:
                errors[i] += (genoProbs[1,i] + genoProbs[2, i] + genoProbs[3, i])
            if genotypes[i] == 1:
                errors[i] += (genoProbs[0, i] + genoProbs[3, i])
            if genotypes[i] == 2:
                errors[i] += (genoProbs[0, i] + genoProbs[1,i] + genoProbs[2, i]) 

def updateSeqError(pedigree, peelingInfo) :
    # The following is a simple EM update for the genotyping error rate at each locus.
    # This update adds the expected number of errors that an individual has marginalizing over their current genotype probabilities.
    # This only uses the homozygotic states, heterozygotic states are ignored (in both the counts + errors terms).
    # We use a max value of 5% and a min value of .0001 percent to make sure the values are reasonable

    counts = np.full(pedigree.nLoci, 1, dtype = np.float32)
    errors = np.full(pedigree.nLoci, 0.001, dtype = np.float32)

    for ind in pedigree:
        if ind.reads is not None:
            updateSeqError_ind(counts, errors, ind.reads[0], ind.reads[1], peelingInfo.getGenoProbs(ind.idn))
    
    newError = errors/counts 
    newError = np.maximum(np.minimum(newError, .01), .0001)
    return newError

@jit(nopython=True)
def updateSeqError_ind(counts, errors, refReads, altReads, genoProbs):
    # Errors occur when genotype is 0 and an alternative allele happens.
    # Errors occur when genotype is 2 (coded as 3) and a reference allele happens.
    # Number of observations is number of reads * probability the individual is homozygous.
    for i in range(len(counts)) :
        counts[i] += (genoProbs[0,i] + genoProbs[3,i]) * (altReads[i] + refReads[i])
        errors[i] += genoProbs[0,i] * altReads[i]
        errors[i] += genoProbs[3,i] * refReads[i]


###
### NOTE: The following code is designed to estimate the recombination rate between markers.
### This is currently broken and does not give good estimates. Strongly recommend not using it, and the associated options have been disabled in tinyPeel.py
###




def updateSeg(peelingInfo):
    #Idea: Split the chromosome into ~10 blocks with slightly different starts/stops. Estimate at each one and then sum.
    #Not sure if this is the best way, but probably will work okay.

    distances = [estDistanceFromBreaks(0, peelingInfo.nLoci-1, nBreaks, peelingInfo) for nBreaks in [5, 10, 20]]
    distance = np.mean(distances[0:2])
    print(distances, ":", distance)

    setupTransmission(distance, peelingInfo)

def estDistanceFromBreaks(loc1, loc2, nBreaks, peelingInfo):
    breakPoints = np.floor(np.linspace(loc1, loc2, nBreaks)).astype(dtype = np.int64)
    distances = [getDistance(breakPoints[i], breakPoints[i+1], peelingInfo) for i in range(nBreaks-1)]
    print(nBreaks, ":", distances)
    return sum(distances)

def getDistance(loc1, loc2, peelingInfo):

    patSeg1 = getSumSeg(loc1, peelingInfo)
    patSeg2 = getSumSeg(loc2, peelingInfo)

    patValid = ((patSeg1 > .99) | (patSeg1 < .01)) & ((patSeg2 > .99) | (patSeg2 < .01)) 
    
    # matValid = ((matSeg1 > .99) | (matSeg1 < .01)) & ((matSeg2 > .99) | (matSeg2 < .01)) 
    
    # patRecomb = np.mean(np.abs(patSeg1[patValid] - patSeg2[patValid]))
    # matRecomb = np.mean(np.abs(matSeg1[matValid] - matSeg2[matValid]))
    
    # recomb = (patRecomb + matRecomb)/2
    
    # difference = np.abs(np.round(patSeg1[patValid]) - np.round(patSeg2[patValid]))
    # return np.mean(difference)

    entropy1 = 1 - (-patSeg1*np.log2(patSeg1) - (1-patSeg1)*np.log2(1-patSeg1))
    entropy2 = 1 - (-patSeg2*np.log2(patSeg2) - (1-patSeg2)*np.log2(1-patSeg2))

    difference = patSeg1*(1-patSeg2) + (1-patSeg1)*patSeg2
    est = np.sum(entropy1*entropy2*difference)/np.sum(entropy1*entropy2)
    # return est
    return haldane(est)

def getSumSeg(loc, peelingInfo):
    seg = peelingInfo.segregation[:, :, loc]
    sumSeg = np.sum(seg, 1)
    patSeg = (seg[:,2] + seg[:,3])/sumSeg
    # matSeg = (seg[:,1] + seg[:,3])/sumSeg
    return patSeg


def haldane(difference) :
    return -np.log(1.0-2.0*difference)/2.0

