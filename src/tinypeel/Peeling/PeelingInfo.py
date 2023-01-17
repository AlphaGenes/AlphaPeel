import concurrent.futures
from numba import jit, float32, int8, int64, optional, boolean
from numba.experimental import jitclass
import numpy as np
from collections import OrderedDict

from ..tinyhouse import InputOutput
from ..tinyhouse import ProbMath
from ..tinyhouse import HaplotypeOperations

import math

#########################################################################
### In this module we define the peeling info object.                 ###  
### This is a just in time container for the various                  ### 
### peeling probability calculations.                                 ### 
#########################################################################

def createPeelingInfo(pedigree, args, createSeg=True, phaseFounder = False) :

    # NOTE: createSeg is added as an option to decrease memory usage during the single locus peeling steps.
    nLoci = pedigree.nLoci

    peelingInfo = jit_peelingInformation(nInd=pedigree.maxIdn, nFam=pedigree.maxFam, nLoci=nLoci, createSeg=createSeg)
    peelingInfo.isSexChrom = args.sexchrom
    # Information about the peeling positions are handled elsewhere.
    peelingInfo.positions = None 

    #Generate the segregation tensors. 
    peelingInfo.segregationTensor = ProbMath.generateSegregation(e = 1e-06)
    peelingInfo.segregationTensor_norm = ProbMath.generateSegregation(e = 1e-06, partial=True) #Partial gives the normalizing constant.

    peelingInfo.segregationTensorXY = ProbMath.generateSegregationXYChrom(e = 1e-06)
    peelingInfo.segregationTensorXX = ProbMath.generateSegregationXXChrom(e = 1e-06)

    peelingInfo.genoError[:] = args.error
    peelingInfo.seqError[:] = args.seqerror
    setupTransmission(args.length, peelingInfo) #Sets up the transmission rates using a custom position list and a total chromosome length.

    for ind in pedigree:
        peelingInfo.sex[ind.idn] = ind.sex

        if ind.genotypes is not None and ind.haplotypes is not None:
            HaplotypeOperations.ind_fillInGenotypesFromPhase(ind)

        sexChromFlag = peelingInfo.isSexChrom and ind.sex == 0 #This is the sex chromosome and the individual is male.

        peelingInfo.penetrance[ind.idn,:,:] = ProbMath.getGenotypeProbabilities(peelingInfo.nLoci, ind.genotypes, ind.reads, peelingInfo.genoError, peelingInfo.seqError, sexChromFlag)
        
        # Set the genotyping/read status for each individual. This will be used for, e.g., estimating the minor allele frequency.
        if ind.genotypes is not None:
            setGenotypeStatusGenotypes(ind.idn, ind.genotypes, peelingInfo)

        if ind.reads is not None:
            setGenotypeStatusReads(ind.idn, ind.reads[0], ind.reads[1], peelingInfo)

        if ind.isGenotypedFounder() and phaseFounder and ind.genotypes is not None:
            loci = getHetMidpoint(ind.genotypes)
            if loci is not None:
                e = args.error
                peelingInfo.penetrance[ind.idn,:,loci] = np.array([e/3, e/3, 1-e, e/3], dtype = np.float32)
    

    if args.penetrance is not None:
        if args.sexchrom:
            print("Using an external penetrance file and the sexchrom option is highly discouraged. Please do not use.")

        if args.esterrors :
            print("External penetrance file included, but esterrors flag used. The two options are incompatible. esterrors set to false.")
            args.esterrors = False
        for pen in args.penetrance:
            addPenetranceFromExternalFile(pedigree, peelingInfo, pen, args)
    # updateMaf(pedigree, peelingInfo)
    return peelingInfo

def setupTransmission(length, peelingInfo) :
    if peelingInfo.positions is None:
        localMap = np.linspace(0, 1, num = peelingInfo.nLoci, dtype = np.float32)
    else:
        localMap = peelingInfo.positions/peelingInfo.positions[-1] #This should be sorted. Need to add in code to check.
    for i in range(peelingInfo.nLoci -1):
        distance = localMap[i+1] - localMap[i]
        distance = distance * length
        peelingInfo.transmissionRate[i] = distance



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
    # This function allows external penetrance files to be read in and added to the gentoype probabilities for an individual.

    print("Reading in penetrance file:", fileName)
    with open(fileName) as f:
        e = 0
        for line in f:
            parts = line.split(); 
            idx = parts[0]; 
            parts = parts[1:]

            if args.startsnp is not None :
                parts = parts[args.startsnp : args.stopsnp+1] #Offset 1 to include stopsnp

            penetranceLine=np.array([float(val) for val in parts], dtype = np.float32)

            if idx not in pedigree.individuals:
                print("Individual", idx, "not found in pedigree. Individual ignored.")
            else:
                ind = pedigree.individuals[idx]
                peelingInfo.penetrance[ind.idn,e,:] *= penetranceLine
                # Normalizing in terms of SNPs seems like a really bad idea.
                # peelingInfo.penetrance[ind.idn,e,:] /= np.sum(peelingInfo.penetrance[ind.idn,e,:], 0) # Normalization added, just in case.
                e = (e+1)%4


@jit(nopython=True)
def getHetMidpoint(geno):
    nLoci = len(geno)
    midpoint = int(nLoci/2)
    index = 0
    changed = False
    while index < nLoci/2:
        if midpoint + index < nLoci:
            if geno[midpoint + index] == 1:
                return midpoint + index
        if midpoint - index >= 0:
            if geno[midpoint - index] == 1:
                return midpoint - index
        index += 1
    return None

spec = OrderedDict()
spec['nInd'] = int64
spec['nFam'] = int64
spec['nLoci'] = int64

spec['isSexChrom'] = boolean
spec['sex'] = int64[:]
spec['genotyped'] = boolean[:,:] #Maybe this should be removed?

# Individual terms: Each will be nInd x 4 x nLoci
spec['anterior'] = float32[:,:,:]
spec['posterior'] = float32[:,:,:]
spec['penetrance'] = float32[:,:,:]
spec['segregation'] = optional(float32[:,:,:])
spec['pointSeg'] = optional(float32[:,:,:]) # I think we don't use this any more. Potentially could be dropped.

# Family terms. Each will be nFam x 4 x nLoci
spec['posteriorSire_minusFam'] = float32[:,:,:]
spec['posteriorDam_minusFam'] = float32[:,:,:]
spec['posteriorSire_new'] = float32[:,:,:]
spec['posteriorDam_new'] = float32[:,:,:]

# Segregation tensors. Each of these will be either 4x4x4x4 or 4x4x4
spec['segregationTensor'] = optional(float32[:,:,:,:])
spec['segregationTensor_norm'] = optional(float32[:,:,:]) # Note: This one is a bit smaller.
spec['segregationTensorXX'] = optional(float32[:,:,:,:])
spec['segregationTensorXY'] = optional(float32[:,:,:,:])

# Marker specific rates:
spec['genoError'] = optional(float32[:])
spec['seqError'] = optional(float32[:])
spec['transmissionRate'] = optional(float32[:])
spec['maf'] = optional(float32[:])

spec['positions'] = optional(float32[:]) # Not sure we use this.
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

        # These are filled in from createPeelingInfo, above.
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

        if createSeg: # Only removes the point seg term since this is not used for single locus peeling.
            self.pointSeg = np.full((self.nInd, 4, self.nLoci), baseValue, dtype = np.float32)
        else:
            self.pointSeg = None

        self.posteriorSire_minusFam = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)
        self.posteriorDam_minusFam = np.full((self.nFam, 4, self.nLoci), baseValue, dtype = np.float32)

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


