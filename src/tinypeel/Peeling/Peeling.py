import concurrent.futures
from numba import jit, float32, int8, int64, optional, boolean
from numba.experimental import jitclass
import numpy as np
from collections import OrderedDict

from ..tinyhouse import InputOutput
from ..tinyhouse import ProbMath
from ..tinyhouse import HaplotypeOperations

import math

# Defining variables for peel up and peel down. Ideally these would be characters, but numba does not support characters.
PEEL_UP = 0
PEEL_DOWN = 1

# This is the main peeling function.
@jit(nopython=True, nogil=True, locals={'e': float32, 'e4':float32, 'e16':float32, 'e1e':float32})
def peel(family, operation, peelingInfo, singleLocusMode) :

    isSexChrom = peelingInfo.isSexChrom

    e = .000001
    e1e = 1-e
    e4 = e/4
    e16 = e/16

    ### Setup local variables from the peeling information container.
    anterior = peelingInfo.anterior
    penetrance = peelingInfo.penetrance
    posterior = peelingInfo.posterior
    segregation = peelingInfo.segregation
    
    pointSeg = peelingInfo.pointSeg
    segregationTensor = peelingInfo.segregationTensor
    segregationTensor_norm = peelingInfo.segregationTensor_norm

    nLoci = peelingInfo.nLoci
    nOffspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn

    #Creating variables here:
    # childToParents: The projection of each child onto the parent genotypes.
    # childSegs: The segregation estimates for a particular child (These are re-used? so need to be stored)
    # allToParents: The projection of each child onto the parental genotypes.
    # parentsMinustChild: The estimate of the parent's genotypes minus the contribution from a specific individual.

    childToParents = np.full((nOffspring, 4, 4, nLoci), 0, dtype = np.float32)
    childSegTensor = np.full((nOffspring, 4, 4, 4, nLoci), 0, dtype = np.float32)
    allToParents = np.full((4, 4, nLoci), 1.0, dtype = np.float32) 
    parentsMinusChild = np.full((nOffspring, 4, 4, nLoci), 1, dtype = np.float32)


    # Some local variables. currentSeg is the segregation estimate of a child (but may be modified).
    currentSeg = np.full((4, nLoci), 1, dtype = np.float32)

    #Construct the joint parent genotypes based on the parent's anterior, penetrance, and posterior terms minus this family.
    
    probSire = anterior[sire,:,:]*penetrance[sire,:,:] * peelingInfo.posteriorSire_minusFam[fam,:,:]
    probDam = anterior[dam,:,:]*penetrance[dam,:,:] * peelingInfo.posteriorDam_minusFam[fam,:,:]

    probSire = probSire/np.sum(probSire, 0)
    probDam = probDam/np.sum(probDam, 0)

    # Einstien sum notation 1: create the joint parental genotypes based on the probabilities for each parent.
    # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)

    jointParents = getJointParents(probSire, probDam)
    jointParents = jointParents/np.sum(np.sum(jointParents, axis = 0), axis = 0)
    jointParents = (1-e)*jointParents + e/16 # There are 4x4 values for each locus in jointparents.

    # There are 4 values for each locus. Normalization is done here so that jointParents is as accurate as possible.
    # We need the posterior terms here for the peeling up step later.
    probSire = probSire*e1e + e4 
    probDam = probDam*e1e + e4 

    # Now construct the parental genotypes based on within-family information. 

    for index in range(nOffspring):
        child = family.offspring[index]
        
        # The child's estimate is the combination of the posterior term and penetrance term for that child.
        # We are estimating the parent's genotypes so the anterior term is ignored to avoid double counting.
        childValues = posterior[child,:,:] * penetrance[child,:,:]
        childValues = childValues/np.sum(childValues, axis = 0)
        childValues = e1e*childValues + e4

        # METHOD 1: Just use the current segregation of the child. 
        currentSeg[:,:] = segregation[child,:,:]
        currentSeg /= np.sum(currentSeg, 0)

        # METHOD 2: Use the segregation estimate of the child minus the contribution at a particular locus.
        # Currently do not recommend using this.
        # if not singleLocusMode :
        #     currentSeg[:,:] = segregation[child,:,:] / pointSeg[child,:,:]
        #     currentSeg /= np.sum(currentSeg, 0)
        # else:
        #     currentSeg[:,:] = segregation[child,:,:]
        
        if isSexChrom and peelingInfo.sex[child] == 0: #0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXY
        if isSexChrom and peelingInfo.sex[child] == 1: #0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXX

        #Einstien sum notation 2: Create the child-specific segregation tensor using the child's currrent segregation estimate.
        # childSegTensor[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
        createChildSegs(segregationTensor, currentSeg, childSegTensor[index,:,:,:,:])
    
        #Einstien sum notation 3: Estimate the parental genotypes based on the child's genotypes and their segregation tensor.
        # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegTensor[index,:,:,:,:], childValues)
        projectChildGenotypes(childSegTensor[index,:,:,:,:], childValues, childToParents[index,:,:,:])
        
    
    # Method 1: estimate the parents genotype and the child-specific posterior terms using iterative normalizing.
    # for i in range(nOffspring) :
    #     parentsMinusChild[i,:,:,:] = jointParents[:,:,:]

    # for i in range(nOffspring):
    #     allToParents *= childToParents[i,:,:,:]
    #     allToParents /= np.sum(np.sum(allToParents, axis = 0), axis=0)

    #     for j in range(nOffspring) :
    #         if i != j :
    #             parentsMinusChild[j,:,:,:] *= childToParents[i,:,:,:]
    #             parentsMinusChild[j,:,:,:] /= np.sum(np.sum(parentsMinusChild[i,:,:,:], axis = 0), axis=0)

    ##
    # Method 2: estimate the parents genotype and the child-specific posterior terms using a log scale.
    ##
    # for i in range(nOffspring) :
    #     parentsMinusChild[i,:,:,:] = np.log(jointParents[:,:,:])
    # allToParents[:,:,:] = 0
    # # #taking out post estimates.
    # for i in range(nOffspring):
    #     log_childToParents = np.log(childToParents[i,:,:,:])
    #     allToParents += log_childToParents
    #     for j in range(nOffspring) :
    #         if i != j :
    #             parentsMinusChild[j,:,:,:] += log_childToParents

    # Method 3: estimate the parents genotype and the child-specific posterior terms using a slightly smarter log scale.

    for i in range(nOffspring) :
        parentsMinusChild[i,:,:,:] = np.log(jointParents[:,:,:])
    allToParents[:,:,:] = 0
    for i in range(nOffspring):
        log_childToParents = np.log(childToParents[i,:,:,:])
        allToParents += log_childToParents
        parentsMinusChild[i,:,:,:] -= log_childToParents # This is done to take away the setimate for an individual child from their parent's posterior term.
    for i in range(nOffspring):
        parentsMinusChild[i,:,:,:] += allToParents

    # Move from a log-scale to a non-log scale and re-normalize.
    allToParents = expNorm2D(allToParents)
    for i in range(nOffspring):
        parentsMinusChild[i,:,:,:] = expNorm2D(parentsMinusChild[i,:,:,:])


    if operation == PEEL_DOWN:
        for i in range(nOffspring):
            child = family.offspring[i]

            #Einstien sum notation 4: Project the parent genotypes down onto the child genotypes.
            # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegTensor[i,:,:,:,:], parentsMinusChild[i,:,:,:])
            projectParentGenotypes(childSegTensor[i,:,:,:,:], parentsMinusChild[i,:,:,:], anterior[child,:,:])
            anterior[child,:,:] /= np.sum(anterior[child,:,:], 0)
    
    if operation == PEEL_UP :
        # Take the allToParents estimate and combine to estimate the sire and dam's posterior estimates (for this family)

        sirePosterior = combineAndReduceAxis1(allToParents, probDam)
        sirePosterior /= np.sum(sirePosterior, axis = 0)
        sirePosterior = sirePosterior*e1e + e4
        peelingInfo.posteriorSire_new[fam,:,:] = sirePosterior

        damPosterior = combineAndReduceAxis0(allToParents, probSire)
        damPosterior /= np.sum(damPosterior, axis = 0)
        damPosterior = damPosterior*e1e + e4
        peelingInfo.posteriorDam_new[fam,:,:] = damPosterior

    if (not singleLocusMode) and (operation == PEEL_DOWN):
        # Estimate the segregation probabilities for each child.

        for i in range(nOffspring):
            # Child values is the same as in the posterior estimation step above.
            child = family.offspring[i]
            childValues = posterior[child,:,:] * penetrance[child,:,:]
            childValues = childValues/np.sum(childValues, axis = 0)
            childValues = e1e*childValues + e4
            
            if isSexChrom and peelingInfo.sex[child] == 0: #0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXY
            if isSexChrom and peelingInfo.sex[child] == 1: #0 for male, 1 for female.
                segregationTensor = peelingInfo.segregationTensorXX

            #Einstien sum notation 5:
            # pointSeg[child,:,:] = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild[i,:,:,:], childValues)
            #Option 1: Estimate without normalizing.
            # estimateSegregation(segregationTensor, parentsMinusChild[i,:,:,:], childValues, pointSeg[child,:,:])
            #Option 2: Estimate with normalizing. I think this is what we want.
            estimateSegregationWithNorm(segregationTensor, segregationTensor_norm, parentsMinusChild[i,:,:,:], childValues, pointSeg[child,:,:])

            segregation[child,:,:] = (1-e)*collapsePointSeg(pointSeg[child,:,:], peelingInfo.transmissionRate) + e/4

#####
##### The following are a large number of "helper" jit functions that replace the einstien sums in the original scripts.
#####


@jit(nopython=True, nogil=True)
def getJointParents(probSire, probDam):
    # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)
    nLoci = probSire.shape[1]
    output = np.full(shape = (4, 4, nLoci), fill_value = 0, dtype = np.float32)
    for a in range(4) :
        for b in range(4) :
            for i in range(nLoci):
                output[a, b, i] = probSire[a,i] * probDam[b,i]
    return output


@jit(nopython=True, nogil=True)
def createChildSegs(segregationTensor, currentSeg, output):
    # childSegs[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
    nLoci = currentSeg.shape[1]
    output[:,:,:,:] = 0
    for a in range(4) :
        for b in range(4) :
            for c in range(4) :
                for d in range(4) :
                    for i in range(nLoci):
                        output[a, b, c, i] += segregationTensor[a, b, c, d]*currentSeg[d,i]
    
    return output


@jit(nopython=True, nogil=True)
def projectChildGenotypes(childSegs, childValues, output):
    # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegs[index,:,:,:,:], childValues)
    nLoci = childSegs.shape[3]
    output[:,:,:] = 0
    for a in range(4) :
        for b in range(4) :
            for c in range(4) :
                for i in range(nLoci):
                    output[a, b, i] += childSegs[a, b, c, i]*childValues[c,i]
    
    return output

@jit(nopython=True, nogil=True)
def projectParentGenotypes(childSegs, parentValues, output):
    # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegs[i,:,:,:,:], parentsMinusChild[i,:,:,:])
    nLoci = childSegs.shape[3]
    output[:,:]=0
    for a in range(4) :
        for b in range(4) :
            for c in range(4) :
                for i in range(nLoci):
                    output[c, i] += childSegs[a,b,c,i]*parentValues[a,b,i]
    
    return output

@jit(nopython=True, nogil=True)
def estimateSegregation(segregationTensor, parentValues, childValues, output):
    # pointSeg[child,:,:] = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild[i,:,:,:], childValues)
    nLoci = childValues.shape[1]
    output[:,:]=0
    for a in range(4) :
        for b in range(4) :
            for c in range(4) :
                for d in range(4) :
                    for i in range(nLoci):
                        output[d, i] += segregationTensor[a,b,c,d]*parentValues[a,b,i]*childValues[c,i]
    return output

@jit(nopython=True, nogil=True)
def estimateSegregationWithNorm(segregationTensor, segregationTensor_norm, parentValues, childValues, output):
    # pointSeg[child,:,:] = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild[i,:,:,:], childValues)
    nLoci = childValues.shape[1]
    output[:,:]=0
    for a in range(4) :
        for b in range(4) :
            for c in range(4) :
                for d in range(4) :
                    for i in range(nLoci):
                        #Check if norm is 0. Otherwise use norm to normalize.
                        if segregationTensor_norm[a, b, c] != 0:
                            output[d, i] += segregationTensor[a,b,c,d]*parentValues[a,b,i]*childValues[c,i]/segregationTensor_norm[a, b, c]
    return output

@jit(nopython=True, nogil=True)
def combineAndReduceAxis1(jointEstimate, parentEstimate):
    # output = np.einsum("abi, bi-> ai", jointEstimate, parentEstimate)
    nLoci = parentEstimate.shape[1]
    output = np.full((4, nLoci), 0, dtype = np.float32)
    for a in range(4):
        for b in range(4) :
            for i in range(nLoci):
                output[a, i] += jointEstimate[a,b,i]*parentEstimate[b, i]
    return output

@jit(nopython=True, nogil=True)
def combineAndReduceAxis0(jointEstimate, parentEstimate):
    # output = np.einsum("abi, ai-> bi", jointEstimate, parentEstimate)
    nLoci = parentEstimate.shape[1]
    output = np.full((4, nLoci), 0, dtype = np.float32)
    for a in range(4):
        for b in range(4) :
            for i in range(nLoci):
                output[b, i] += jointEstimate[a,b,i]*parentEstimate[a, i]
    return output

@jit(nopython=True, nogil=True)
def expNorm2D(mat):
    # Matrix is 4x4xnLoci: Output is to take the exponential of the matrix and normalize each locus. We need to make sure that there are not any overflow values.
    nLoci = mat.shape[2]
    for i in range(nLoci):
        maxVal = 1 # Log of anything between 0-1 will be less than 0. Using 1 as a default.
        for a in range(4):
            for b in range(4):
                if mat[a, b, i] > maxVal or maxVal == 1:
                    maxVal = mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] -= maxVal
    # Normalize.
    tmp = np.exp(mat)
    for i in range(nLoci):
        total = 0
        for a in range(4):
            for b in range(4):
                total += tmp[a, b, i]
        for a in range(4):
            for b in range(4):
                tmp[a, b, i] /= total
    return tmp



@jit(nopython=True, nogil=True)
def expNorm1D(mat):
    # Matrix is 4x4xnLoci: Output is to take the exponential of the matrix and normalize each locus. We need to make sure that there are not any overflow values.
    nLoci = mat.shape[1]
    for i in range(nLoci):
        maxVal = 1 # Log of anything between 0-1 will be less than 0. Using 1 as a default.
        for a in range(4):
            if mat[a, i] > maxVal or maxVal == 1:
                maxVal = mat[a, i]
        for a in range(4):
            mat[a, i] -= maxVal
    tmp = np.exp(mat)
    for i in range(nLoci):
        total = 0
        for a in range(4):
            total += tmp[a,i]
        for a in range(4):
            tmp[a, i] /= total
    return tmp


@jit(nopython=True, nogil=True, locals={'e': float32, 'e2':float32, 'e1e':float32, 'e2i':float32})
def collapsePointSeg(pointSeg, transmission):

    # This is the forward backward algorithm.
    # Segregation estimate state ordering: pp, pm, mp, mm
    nLoci = pointSeg.shape[1] 

    seg = np.full(pointSeg.shape, .25, dtype = np.float32)
    for i in range(nLoci):
        for j in range(4):
            seg[j,i] = pointSeg[j,i]

    tmp = np.full((4), 0, dtype = np.float32)
    new = np.full((4), 0, dtype = np.float32)

    prev = np.full((4), .25, dtype = np.float32)
    for i in range(1, nLoci):
        e = transmission[i-1]
        e2 = e**2
        e1e = e*(1-e)
        e2i = (1.0-e)**2
        for j in range(4):
            tmp[j] = prev[j]*pointSeg[j,i-1]
        
        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j]/sum_j

        # !                  fm  fm  fm  fm 
        # !segregationOrder: pp, pm, mp, mm

        new[0] = e2*tmp[3] + e1e*(tmp[1] + tmp[2]) + e2i*tmp[0] 
        new[1] = e2*tmp[2] + e1e*(tmp[0] + tmp[3]) + e2i*tmp[1] 
        new[2] = e2*tmp[1] + e1e*(tmp[0] + tmp[3]) + e2i*tmp[2] 
        new[3] = e2*tmp[0] + e1e*(tmp[1] + tmp[2]) + e2i*tmp[3] 

        # tmp = tmp/np.sum(tmp)
        # new = e2i*tmp + e2 + e1e*(tmp[0] + tmp[3])*same + e1e*(tmp[1] + tmp[2])*diff       

        for j in range(4):
            seg[j,i] *= new[j]
        # seg[:,i] *= new
        prev = new

    prev = np.full((4), .25, dtype = np.float32)
    for i in range(nLoci-2, -1, -1): #zero indexed then minus one since we skip the boundary.
        e = transmission[i]
        e2 = e**2
        e1e = e*(1-e)
        e2i = (1.0-e)**2
        
        for j in range(4):
            tmp[j] = prev[j]*pointSeg[j,i+1]
        
        sum_j = 0
        for j in range(4):
            sum_j += tmp[j]
        for j in range(4):
            tmp[j] = tmp[j]/sum_j

        new[0] = e2*tmp[3] + e1e*(tmp[1] + tmp[2]) + e2i*tmp[0] 
        new[1] = e2*tmp[2] + e1e*(tmp[0] + tmp[3]) + e2i*tmp[1] 
        new[2] = e2*tmp[1] + e1e*(tmp[0] + tmp[3]) + e2i*tmp[2] 
        new[3] = e2*tmp[0] + e1e*(tmp[1] + tmp[2]) + e2i*tmp[3] 

        for j in range(4):
            seg[j,i] *= new[j]
        prev = new
    
    for i in range(nLoci):
        sum_j = 0
        for j in range(4):
            sum_j += seg[j, i]
        for j in range(4):
            seg[j, i] = seg[j, i]/sum_j

    # seg = seg/np.sum(seg, 0)
    return(seg)
