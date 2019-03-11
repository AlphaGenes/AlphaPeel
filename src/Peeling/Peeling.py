import concurrent.futures
from numba import jit, jitclass, float32, int8, int64, optional, boolean

import numpy as np
from collections import OrderedDict
from tinyhouse import InputOutput
from tinyhouse import ProbMath

import math

# @jit(nopython=True)
PEEL_UP = 0
PEEL_DOWN = 1

@jit(nopython=True, nogil=True, locals={'e': float32, 'e4':float32, 'e16':float32, 'e1e':float32, 'childValues':float32[:,:]})
def peel(family, operation, peelingInfo, singleLocusMode) :

    isSexChrom = peelingInfo.isSexChrom

    e = .000001
    e1e = 1-e
    e4 = e/4
    e16 = e/16
    ### Setup
    anterior = peelingInfo.anterior
    penetrance = peelingInfo.penetrance
    posterior = peelingInfo.posterior
    segregation = peelingInfo.segregation
    
    pointSeg= peelingInfo.pointSeg
    segregationTensor = peelingInfo.segregationTensor
    segregationTensor_norm = peelingInfo.segregationTensor_norm

    nLoci = peelingInfo.nLoci
    nOffspring = len(family.offspring)
    sire = family.sire
    dam = family.dam
    fam = family.idn
    #Creating variables here:


    childToParents = np.full((nOffspring, 4, 4, nLoci), 0, dtype = np.float32)
    childSegs = np.full((nOffspring, 4, 4, 4, nLoci), 0, dtype = np.float32)
    allToParents = np.full((4, 4, nLoci), 1.0, dtype = np.float32) 
    parentsMinusChild = np.full((nOffspring, 4, 4, nLoci), 1, dtype = np.float32)

    currentSeg = np.full((4, nLoci), 1, dtype = np.float32)
    nullSeg = np.full(currentSeg.shape, .25, dtype = np.float32)
    nullSegTensor = np.full((4, 4, 4, nLoci), .25, dtype = np.float32)

    #Actual code.
    

    # probSire = anterior[sire,:,:]*penetrance[sire,:,:]*posterior[sire,:,:] / peelingInfo.posteriorSire[fam,:,:]
    # probDam = anterior[dam,:,:]*penetrance[dam,:,:]*posterior[dam,:,:] / peelingInfo.posteriorDam[fam,:,:]

    probSire = anterior[sire,:,:]*penetrance[sire,:,:] * peelingInfo.posteriorSire[fam,:,:]
    probDam = anterior[dam,:,:]*penetrance[dam,:,:] * peelingInfo.posteriorDam[fam,:,:]

    probSire = probSire/np.sum(probSire, 0)
    probDam = probDam/np.sum(probDam, 0)

    #Einstien sum notation 1:
    # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)
    # jointParents = probSire[:,None,:]*probDam[None,:,:]
    jointParents = getJointParents(probSire, probDam)
    jointParents = jointParents/np.sum(np.sum(jointParents, axis = 0), axis = 0)
    jointParents = (1-e)*jointParents + e/16 #Each loci has 16 values.

    probSire = probSire*e1e + e4
    probDam = probDam*e1e + e4

    createChildSegs(segregationTensor, nullSeg, nullSegTensor)

    for index in range(nOffspring):
        child = family.offspring[index]
        
        childValues = posterior[child,:,:] * penetrance[child,:,:]
        childValues = childValues/np.sum(childValues, axis = 0)
        childValues = e1e*childValues + e4

        #METHOD 1    
        # if not singleLocusMode :
        #     currentSeg[:,:] = segregation[child,:,:] / pointSeg[child,:,:]
        #     currentSeg /= np.sum(currentSeg, 0)
        # else:
        #     currentSeg[:,:] = segregation[child,:,:]
        #METHOD 2
        currentSeg[:,:] = segregation[child,:,:]
        currentSeg /= np.sum(currentSeg, 0)
        
        # currentSeg = nullSeg
        if isSexChrom and peelingInfo.sex[child] == 0: #0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXY
        if isSexChrom and peelingInfo.sex[child] == 1: #0 for male, 1 for female.
            segregationTensor = peelingInfo.segregationTensorXX

        #Einstien sum notation 2:
        # childSegs[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
        # childSegs[index,:,:,:,:] = np.sum(currentSeg[None,None,None,:,:]*segregationTensor[:,:,:,:,None], 3)
        createChildSegs(segregationTensor, currentSeg, childSegs[index,:,:,:,:])
    
        #Einstien sum notation 3:
        # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegs[index,:,:,:,:], childValues)
        # childToParents[index,:,:,:] = np.sum(childSegs[index,:,:,:,:] * childValues[None,None,:,:], 2)
        projectChildGenotypes(childSegs[index,:,:,:,:], childValues, childToParents[index,:,:,:])
        
        # disable peeling up for going to parents.
        # projectChildGenotypes(nullSegTensor, childValues, childToParents[index,:,:,:])
        # print(childToParents[index,:,:,:])


    #Now collapse to get the current posterior estimate for the parents inc. all children.
    
    # #Method 1: use iterative normalizing.
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
    #Method 2: Log scale + final conversion.
    ##
    for i in range(nOffspring) :
        parentsMinusChild[i,:,:,:] = np.log(jointParents[:,:,:])
    allToParents[:,:,:] = 0
    # #taking out post estimates.
    for i in range(nOffspring):
        log_childToParents = np.log(childToParents[i,:,:,:])
        allToParents += log_childToParents
        for j in range(nOffspring) :
            if i != j :
                parentsMinusChild[j,:,:,:] += log_childToParents
    
    #This drops everyone to be back on the right scale, and normalizes.
    allToParents = expNorm(allToParents)
    allToParents /= np.sum(np.sum(allToParents, axis = 0), axis=0)
    for i in range(nOffspring):
        parentsMinusChild[i,:,:,:] = expNorm(parentsMinusChild[i,:,:,:])
        parentsMinusChild[i,:,:,:] /= np.sum(np.sum(parentsMinusChild[i,:,:,:], axis = 0), axis=0)


    # print(parentsMinusChild[0,:,:,0])
    ###


    #Now perform operations. 
    if operation == PEEL_DOWN:
        for i in range(nOffspring):
            child = family.offspring[i]
            #Einstien sum notation 4:
            # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegs[i,:,:,:,:], parentsMinusChild[i,:,:,:])
            # anterior[child,:,:] = np.sum(childSegs[i,:,:,:,:]*parentsMinusChild[i,:,:,None,:], (0,1))
            projectParentGenotypes(childSegs[i,:,:,:,:], parentsMinusChild[i,:,:,:], anterior[child,:,:])
            # print(anterior[child,:,:] )

            anterior[child,:,:] /= np.sum(anterior[child,:,:], 0)
    
    if operation == PEEL_UP :
        # sirePosterior = np.sum(allToParents, axis =1)*(1-e) + e/4
        sirePosterior = combineAndReduceAxis1(allToParents, probDam)
        sirePosterior /= np.sum(sirePosterior, axis = 0)
        sirePosterior = sirePosterior*e1e + e4
        peelingInfo.posteriorSire_new[fam,:,:] = sirePosterior

        # damPosterior = np.sum(allToParents, axis =0)*(1-e) + e/4
        damPosterior = combineAndReduceAxis0(allToParents, probSire)
        damPosterior /= np.sum(damPosterior, axis = 0)
        damPosterior = damPosterior*e1e + e4

        peelingInfo.posteriorDam_new[fam,:,:] = damPosterior

    if (not singleLocusMode) and (operation == PEEL_DOWN):
        for i in range(nOffspring):
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
            # tmp = np.sum(segregationTensor[:,:,:,:,None] * parentsMinusChild[i,:,:,None,None,:],(0,1))
            # pointSeg[child,:,:] = np.sum(tmp[:,:,:] * childValues[:,None,:], 0)
            #Option 1
            # estimateSegregation(segregationTensor, parentsMinusChild[i,:,:,:], childValues, pointSeg[child,:,:])
            #Option 2
            estimateSegregationWithNorm(segregationTensor, segregationTensor_norm, parentsMinusChild[i,:,:,:], childValues, pointSeg[child,:,:])

            # print(pointSeg[child,:,:] )
            # segregation[child,:,:] = (1-e)*collapsePointSeg2(pointSeg[child,:,:], peelingInfo.transmissionRate) + e/4
            segregation[child,:,:] = (1-e)*collapsePointSeg3(pointSeg[child,:,:], peelingInfo.transmissionRate) + e/4

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
    nLoci = parentEstimate.shape[1]
    output = np.full((4, nLoci), 0, dtype = np.float32)
    for a in range(4):
        for b in range(4) :
            for i in range(nLoci):
                output[a, i] += jointEstimate[a,b,i]*parentEstimate[b, i]
    return output

@jit(nopython=True, nogil=True)
def combineAndReduceAxis0(jointEstimate, parentEstimate):
    nLoci = parentEstimate.shape[1]
    output = np.full((4, nLoci), 0, dtype = np.float32)
    for a in range(4):
        for b in range(4) :
            for i in range(nLoci):
                output[b, i] += jointEstimate[a,b,i]*parentEstimate[a, i]
    return output

@jit(nopython=True, nogil=True)
def expNorm(mat):
    # Matrix is 4x4xnLoci
    nLoci = mat.shape[2]
    for i in range(nLoci):
        maxVal = 1 #Log of anything between 0-1 will be less than 0. Using 1 as a default.
        for a in range(4):
            for b in range(4):
                if mat[a, b, i] > maxVal or maxVal == 1:
                    maxVal = mat[a, b, i]
        for a in range(4):
            for b in range(4):
                mat[a, b, i] -= maxVal
    return np.exp(mat)

@jit(nopython=True, nogil=True)
def expNorm1(mat):
    # Matrix is 4xnLoci
    nLoci = mat.shape[1]
    for i in range(nLoci):
        maxVal = 1
        for a in range(4):
            if mat[a, i] > maxVal or maxVal == 1:
                maxVal = mat[a, i]
        for a in range(4):
            mat[a, i] -= maxVal
    return np.exp(mat)


##########
##########
##########
##########

# @jit(nopython=True, nogil=True)
# def getJointParents(probSire, probDam):
#     # jointParents = np.einsum("ai, bi -> abi", probSire, probDam)
#     nLoci = probSire.shape[1]
#     output = np.full(shape = (4, 4, nLoci), fill_value = 0, dtype = np.float32)
#     for a in range(4) :
#         for b in range(4) :
#             output[a, b, :] = probSire[a,:] * probDam[b,:]
#     return output


# @jit(nopython=True, nogil=True)
# def createChileSegs(segregationTensor, currentSeg, output):
#     # childSegs[index,:,:,:,:] = np.einsum("abcd, di -> abci", segregationTensor, currentSeg)
#     output[:,:,:,:] = 0
#     for a in range(4) :
#         for b in range(4) :
#             for c in range(4) :
#                 for d in range(4) :
#                     output[a, b, c, :] += segregationTensor[a, b, c, d]*currentSeg[d,:]
    
#     return output


# @jit(nopython=True, nogil=True)
# def projectChildGenotypes(childSegs, childValues, output):
#     # childToParents[index,:,:,:] = np.einsum("abci, ci -> abi", childSegs[index,:,:,:,:], childValues)
#     output[:,:,:] = 0
#     for a in range(4) :
#         for b in range(4) :
#             for c in range(4) :
#                 output[a, b, :] += childSegs[a, b, c, :]*childValues[c,:]
    
#     return output

# @jit(nopython=True, nogil=True)
# def projectParentGenotypes(childSegs, parentValues, output):
#     # anterior[child,:,:] = np.einsum("abci, abi -> ci", childSegs[i,:,:,:,:], parentsMinusChild[i,:,:,:])
#     output[:,:]=0
#     for a in range(4) :
#         for b in range(4) :
#             for c in range(4) :
#                 output[c, :] += childSegs[a,b,c,:]*parentValues[a,b,:]
    
#     return output

# @jit(nopython=True, nogil=True)
# def estimateSegregation(segregationTensor, parentValues, childValues, output):
#     # pointSeg[child,:,:] = np.einsum("abcd, abi, ci-> di", segregationTensor, parentsMinusChild[i,:,:,:], childValues)
#     output[:,:]=0
#     for a in range(4) :
#         for b in range(4) :
#             for c in range(4) :
#                 for d in range(4) :
#                     output[d, :] += segregationTensor[a,b,c,d]*parentValues[a,b,:]*childValues[c,:]
#     return output


##########
##########
##########
##########


# @jit(nopython=True)
# def collapsePointSeg(pointSeg, transmission):
#     nLoci = pointSeg.shape[1] 
#     workLeft = np.full(pointSeg.shape, .25, dtype = np.float32)
#     workRight = np.full(pointSeg.shape, .25, dtype = np.float32)

#     ##Note: We can probably be faaaaaaar more clever here and not do the matrix mult.
#     for i in range(1, nLoci):
#         workLeft[:,i] = workLeft[:,i-1]*pointSeg[:,i-1]
#         workLeft[:,i] = np.dot(transmission[:,:, i-1], workLeft[:,i]/np.sum(workLeft[:,i]))
#     for i in range(nLoci-2, -1, -1): #zero indexed then minus one since we skip the boundary.
#         workRight[:,i] = workRight[:,i+1]*pointSeg[:,i+1]
#         workRight[:,i] = np.dot(transmission[:,:,i], workRight[:,i]/np.sum(workRight[:,i]))

#     seg = workLeft*workRight*pointSeg

#     seg = seg/np.sum(seg, 0)
#     return(seg)

@jit(nopython=True, nogil=True, locals={'e': float32, 'e2':float32, 'e1e':float32, 'e2i':float32})
def collapsePointSeg2(pointSeg, transmission):
    nLoci = pointSeg.shape[1] 

    seg = np.full(pointSeg.shape, .25, dtype = np.float32)
    seg[:,:] = pointSeg[:,:]

    same = np.array((0,1,1,0), dtype = np.float32)
    diff = np.array((1,0,0,1), dtype = np.float32)
    tmp = np.full((4), 0, dtype = np.float32)
    new = np.full((4), 0, dtype = np.float32)

    prev = np.full((4), .25, dtype = np.float32)
    for i in range(1, nLoci):
        e = transmission[i-1]
        e2 = e**2
        e1e = e*(1-e)
        e2i = 1.0 - e2

        tmp = prev*pointSeg[:,i-1]
        tmp = tmp/np.sum(tmp)
        new = e2i*tmp + e2 + e1e*(tmp[0] + tmp[3])*same + e1e*(tmp[1] + tmp[2])*diff       
        seg[:,i] *= new
        prev = new

    prev = np.full((4), .25, dtype = np.float32)
    for i in range(nLoci-2, -1, -1): #zero indexed then minus one since we skip the boundary.
        e = transmission[i]
        e2 = e**2
        e1e = e*(1-e)
        e2i = 1.0 - e2

        tmp = prev*pointSeg[:,i+1]
        tmp = tmp/np.sum(tmp)

        new = e2i*tmp + e2 + e1e*(tmp[0] + tmp[3])*same + e1e*(tmp[1] + tmp[2])*diff       
        seg[:,i] *= new
        prev = new

    seg = seg/np.sum(seg, 0)
    return(seg)

@jit(nopython=True, nogil=True, locals={'e': float32, 'e2':float32, 'e1e':float32, 'e2i':float32})
def collapsePointSeg3(pointSeg, transmission):
    #pp, pm, mp, mm
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

### END ACTUAL MATH

###UPDATE TERMS
# @profile
def updateMaf(pedigree, peelingInfo):
    print("Estimating maf")
    maf = np.full(peelingInfo.nLoci, .5, dtype = np.float32)
    for i in range(peelingInfo.nLoci):
        # print("Maf update", i)
        # mafRegion = np.full(nBreaks - 1, 0, dtype = np.float32)
        # for j in range(nBreaks -1):
        #     mafRegion[j] = mafLoglikelihood(peelingInfo, 1.0/nBreaks*(j+1), i)
        # maf[i] = 1.0/nBreaks*(np.argmax(mafRegion) + 1)
        maf[i] = newtonMafUpdates(peelingInfo, i)

    mafGeno = ProbMath.getGenotypesFromMaf(maf)
    # print(mafGeno)
    for ind in pedigree:
        if ind.isFounder():
            peelingInfo.anterior[ind.idn,:,:] = mafGeno
    peelingInfo.maf = maf.astype(np.float32)

def newtonMafUpdates(peelingInfo, index) :
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
    #Log likelihood first + second derivitives
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

@jit(nopython = True)
def mafLoglikelihood(peelingInfo, maf, index):
    score = 0
    maf_squared = maf**2
    maf_one_minus_maf = maf*(1-maf)
    one_minus_maf_squared = (1-maf)**2

    for i in range(peelingInfo.nInd):
        if peelingInfo.genotyped[i, index] :
        # if True :
            genoProbs = peelingInfo.penetrance[i,:,index]

            prob = 0
            prob += one_minus_maf_squared*genoProbs[0]
            prob += maf_one_minus_maf*genoProbs[1]
            prob += maf_one_minus_maf*genoProbs[2]
            prob += maf_squared*genoProbs[3]
            score += math.log(prob)
    return score



def updatePenetrance(pedigree, peelingInfo):
    peelingInfo.genoError = updateGenoError(pedigree, peelingInfo)
    peelingInfo.seqError = updateSeqError(pedigree, peelingInfo)

    if peelingInfo.isSexChrom:
        print("Updating error rates and minor allele frequencies for sex chromosomes are not well test and will break in interesting ways.")
    for ind in pedigree:
        sexChromFlag = peelingInfo.isSexChrom and ind.sex == 0 #This is the sex chromosome and the individual is male.
        peelingInfo.penetrance[ind.idn,:,:] = ProbMath.getGenotypeProbabilities(peelingInfo.nLoci, ind.genotypes, ind.reads, peelingInfo.genoError, peelingInfo.seqError, sexChromFlag)


def updateGenoError(pedigree, peelingInfo) :
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
        if genotypes[i] != 9:
            counts[i] += 1
            if genotypes[i] == 0:
                errors[i] += (genoProbs[1,i] + genoProbs[2, i] + genoProbs[3, i])
            if genotypes[i] == 1:
                errors[i] += (genoProbs[0, i] + genoProbs[3, i])
            if genotypes[i] == 2:
                errors[i] += (genoProbs[0, i] + genoProbs[1,i] + genoProbs[2, i]) 

def updateSeqError(pedigree, peelingInfo) :
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
    for i in range(len(counts)) :
        counts[i] += (genoProbs[0,i] + genoProbs[3,i]) * (altReads[i] + refReads[i])
        errors[i] += genoProbs[0,i] * altReads[i]
        errors[i] += genoProbs[3,i] * refReads[i]

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

def setupTransmission(length, peelingInfo) :
    if peelingInfo.positions is None:
        localMap = np.linspace(0, 1, num = peelingInfo.nLoci, dtype = np.float32)
    else:
        localMap = peelingInfo.positions/peelingInfo.positions[-1] #This should be sorted. Need to add in code to check.
    for i in range(peelingInfo.nLoci -1):
        distance = localMap[i+1] - localMap[i]
        distance = distance * length
        peelingInfo.transmissionRate[i] = distance

###PEELING OBJECT CREATION

def createPeelingInfo(pedigree, args, createSeg=True, start = None, stop = None, phaseFounder = False) :
    if start is None: start = 0
    if stop is None: stop = pedigree.nLoci
    nLoci = stop - start
    print(start, stop)
    peelingInfo = jit_peelingInformation(nInd=pedigree.maxIdn, nFam=pedigree.maxFam, nLoci=nLoci, createSeg=createSeg)
    
    peelingInfo.isSexChrom = args.sexchrom

    peelingInfo.positions = None 
    # peelingInfo.positions = InputOutput.getPositionsFromMap(args.mapfile, peelingInfo.nLoci)

    #Generate the segregation matrix
    peelingInfo.segregationTensor = ProbMath.generateSegregation(e = 1e-06)
    peelingInfo.segregationTensor_norm = ProbMath.generateSegregation(e = 1e-06, partial=True) #Partial gives the normalizing constant.

    peelingInfo.segregationTensorXY = ProbMath.generateSegregationXYChrom(e = 1e-06)
    peelingInfo.segregationTensorXX = ProbMath.generateSegregationXXChrom(e = 1e-06)

    peelingInfo.genoError[:] = args.error
    peelingInfo.seqError[:] = args.seqerror

    #setup Transmission Probabilities, error rates, maf

    setupTransmission(args.length, peelingInfo)
    count = 0



    for ind in pedigree:
        peelingInfo.sex[ind.idn] = ind.sex

        # if ind.genotypes is not None and ind.haplotypes is not None:
        #     Imputation.ind_fillInGenotypesFromPhase(ind)

        if ind.genotypes is None: genotypes = None
        else: genotypes = ind.genotypes[start:stop]

        if ind.reads is None: reads = None
        else: reads = (ind.reads[0][start:stop], ind.reads[1][start:stop])

        sexChromFlag = peelingInfo.isSexChrom and ind.sex == 0 #This is the sex chromosome and the individual is male.
        if genotypes is not None:
            setGenotypeStatusGenotypes(ind.idn, genotypes, peelingInfo)


        if reads is not None:
            setGenotypeStatusReads(ind.idn, reads[0], reads[1], peelingInfo)

        peelingInfo.penetrance[ind.idn,:,:] = ProbMath.getGenotypeProbabilities(peelingInfo.nLoci, genotypes, reads, peelingInfo.genoError, peelingInfo.seqError, sexChromFlag)

        

        if ind.isGenotypedFounder() and phaseFounder and genotypes is not None:
            loci = getHetMidpoint(genotypes)
            count = count + 1
            # print(count, ind.idn, ind.idx, loci)

            if loci is not None:
                e = args.error
                peelingInfo.penetrance[ind.idn,:,loci] = np.array([e/3, e/3, 1-e, e/3], dtype = np.float32)
    

    if args.penetrance is not None:
        if args.sexchrom:
            print("Using an external penetrance file and the sexchrom option is highly discouraged. Please do not use unless absolutely certain.")

        if args.esterrors :
            print("External penetrance file included, but esterrors flag used. The two options are incompatible. esterrors set to false.")
            args.esterrors = False
        for pen in args.penetrance:
            addPenetranceFromExternalFile(pedigree, peelingInfo, pen, args)
    # updateMaf(pedigree, peelingInfo)
    return peelingInfo

@jit(nopython=True)
def setGenotypeStatusGenotypes(idn, genotypes, peelingInfo):
    nLoci = len(genotypes)
    if genotypes is not None:
        for i in range(nLoci):
            peelingInfo.genotyped[idn, i] = peelingInfo.genotyped[idn, i] or genotypes[i] != 9

@jit(nopython=True)
def setGenotypeStatusReads(idn, reads0, reads1, peelingInfo):
    nLoci = len(reads0)
    if reads0 is not None and reads1 is not None:
        for i in range(nLoci):
            peelingInfo.genotyped[idn, i] = peelingInfo.genotyped[idn, i] or reads0[i] != 0 or reads1[i] != 0

def addPenetranceFromExternalFile(pedigree, peelingInfo, fileName, args):

    print("Reading in penetrance file:", fileName)
    
    with open(fileName) as f:
        e = 0
        for line in f:
            parts = line.split(); 
            idx = parts[0]; 
            parts = parts[1:]

            if args.startsnp is not None :
                parts = parts[args.startsnp : args.stopsnp + 1] #Offset 1 for id and 2 for id + include stopsnp

            penetranceLine=np.array([float(val) for val in parts], dtype = np.float32)

            if idx not in pedigree.individuals:
                print("Individual", idx, "not found in pedigree. Individual ignored.")
            else:
                ind = pedigree.individuals[idx]
                peelingInfo.penetrance[ind.idn,e,:] *= penetranceLine
                e = (e+1)%4


@jit(nopython=True)
def getHetMidpoint(geno):
    midpoint = int(len(geno)/2)
    index = 0
    e = 1
    changed = False
    while not changed :
        if geno[midpoint + index * e] == 1:
            return midpoint + index * e
        e = -e
        if e == -1: index += 1
        if index >= midpoint: changed = True
    return None

spec = OrderedDict()
spec['nInd'] = int64
spec['nFam'] = int64
spec['nLoci'] = int64

spec['isSexChrom'] = boolean
spec['sex'] = int64[:]
spec['genotyped'] = boolean[:,:]

spec['anterior'] = float32[:,:,:]
spec['posterior'] = float32[:,:,:]
spec['penetrance'] = float32[:,:,:]
spec['segregation'] = optional(float32[:,:,:]) # nInd, 4, nLoci
spec['pointSeg'] = optional(float32[:,:,:])

spec['posteriorSire'] = float32[:,:,:]
spec['posteriorDam'] = float32[:,:,:]
spec['posteriorSire_new'] = float32[:,:,:]
spec['posteriorDam_new'] = float32[:,:,:]

# spec['epsilon'] = optional(float32[:,:,:])
spec['segregationTensor'] = optional(float32[:,:,:,:])
spec['segregationTensorXX'] = optional(float32[:,:,:,:])
spec['segregationTensorXY'] = optional(float32[:,:,:,:])

spec['segregationTensor_norm'] = optional(float32[:,:,:])


spec['genoError'] = optional(float32[:])
spec['seqError'] = optional(float32[:])
spec['transmissionRate'] = optional(float32[:])
spec['maf'] = optional(float32[:])

spec['positions'] = optional(float32[:])

spec['iteration'] = int64
@jitclass(spec)
class jit_peelingInformation(object):
    def __init__(self, nInd, nFam, nLoci, createSeg=True):
        self.iteration = 0
        self.nInd = nInd
        self.nFam = nFam
        self.nLoci = nLoci
        
        self.isSexChrom = False

        self.construct(createSeg)
        # self.epsilon = None
        self.positions = None
        self.segregationTensor = None
        self.segregationTensor_norm = None

    def construct(self, createSeg = True) :
        baseValue = .25
        self.sex = np.full(self.nInd, 0, dtype = np.int64)

        self.genotyped = np.full((self.nInd, self.nLoci), False, dtype = np.bool_)

        self.anterior = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)
        self.posterior = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)
        self.penetrance = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)
        
        self.segregation = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)

        if createSeg:
            self.pointSeg = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)
        else:
            self.pointSeg = None

        self.posteriorSire = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)
        self.posteriorDam = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)

        self.posteriorSire_new = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)
        self.posteriorDam_new = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)

        self.genoError = np.full((self.nLoci), 0, dtype = np.float32)
        self.seqError = np.full((self.nLoci), 0, dtype = np.float32)
        self.maf = np.full((self.nLoci), .5, dtype = np.float32)
        self.transmissionRate = np.full((self.nLoci-1), 0, dtype = np.float32)


    def getGenoProbs(self, idn):
        genoProbs = self.anterior[idn,:,:]*self.posterior[idn,:,:]*self.penetrance[idn,:,:]
        genoProbs = genoProbs/np.sum(genoProbs,0) 
        return genoProbs

    # def getDosages(self, idn):
    #     genoProbs = self.getGenoProbs(idn)
    #     value = np.dot(np.array([0,1,1,2], dtype = np.float32), genoProbs)

    #     return 
    # def getDosagesForMaf(self, idn):
    #     genoProbs = self.posterior[idn,:,:]*self.penetrance[idn,:,:]
    #     # genoProbs = self.penetrance[idn,:,:]
    #     genoProbs = genoProbs/np.sum(genoProbs,0) 
    #     return np.dot(np.array([0,1,1,2], dtype = np.float32), genoProbs)




#I'm 90% sure we don't need to do this... Oh. The transmission structure on this is more complicated.
#Yeah, we need this. Let me sort this out later.

