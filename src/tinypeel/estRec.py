import numpy as np
from numba import jit, jitclass, float32, int32, int64, optional

from .tinyhouse import Pedigree
from .tinyhouse import ProbMath
from .tinyhouse import InputOutput 
from .tinyhouse import HaplotypeOperations 

from .Peeling import Peeling
from .Peeling import PeelingIO
from .Peeling import PeelingInfo
from .Peeling import PeelingUpdates

import concurrent.futures
from itertools import repeat



def estimateRecombinations(pedigree, peelingInfo, args):
    nWorkers = args.maxthreads
    singleLocusMode = False
    no_post = True

    for families in reversed(pedigree.genFamilies):
        print("Peeling Up")
        jit_families = [family.toJit() for family in families]

        for family in families:
            fill(family.sire, peelingInfo, phase=True)
            fill(family.dam, peelingInfo, phase=True)
            for child in family.offspring:
                fill(child, peelingInfo, phase=False)


        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers) as executor:
                 results = executor.map(Peeling.peel, jit_families, repeat(Peeling.PEEL_UP), repeat(peelingInfo), repeat(singleLocusMode), repeat(no_post))
        else:
            for family in jit_families:
                Peeling.peel(family, Peeling.PEEL_UP, peelingInfo, singleLocusMode, args.no_post)

    mapLength = Peeling.setRecombinations(pedigree, peelingInfo)

def fill(ind, peelingInfo, phase):
    if phase:
        updatePenetrance_phase(ind.haplotypes[0], ind.haplotypes[1], peelingInfo.penetrance[ind.idn,:,:])
    else:
        updatePenetrance_no_phase(ind.haplotypes[0], ind.haplotypes[1], peelingInfo.penetrance[ind.idn,:,:])

def updatePenetrance_phase(hap0, hap1, penetrance):
    nLoci = len(hap0)
    e = 0.01
    for i in range(nLoci):
        penetrance[:,i] = [.25, .25, .25, .25]
        if hap0[i] != 9 and hap1[i] != 9:
            if hap0[i] == 0 and hap1[i] == 0:
                penetrance[:,i] = [1-e, e/3, e/3, e/3]
            
            if hap0[i] == 0 and hap1[i] == 1:
                penetrance[:,i] = [e/3, 1-e, e/3, e/3]
            
            if hap0[i] == 1 and hap1[i] == 0:
                penetrance[:,i] = [e/3, e/3, 1-e, e/3]
            
            if hap0[i] == 1 and hap1[i] == 1:
                penetrance[:,i] = [e/3, e/3, e/3, 1-e]

def updatePenetrance_no_phase(hap0, hap1, penetrance):
    nLoci = len(hap0)
    e = 0.01
    for i in range(nLoci):
        penetrance[:,i] = [.25, .25, .25, .25]
        if hap0[i] != 9 and hap1[i] != 9:
            if hap0[i] == 0 and hap1[i] == 0:
                penetrance[:,i] = [1-e, e/3, e/3, e/3]
            
            if hap0[i] == 0 and hap1[i] == 1:
                penetrance[:,i] = [e/2, (1-e)/2, (1-e)/2, e/2]
            
            if hap0[i] == 1 and hap1[i] == 0:
                penetrance[:,i] = [e/2, (1-e)/2, (1-e)/2, e/2]
            
            if hap0[i] == 1 and hap1[i] == 1:
                penetrance[:,i] = [e/3, e/3, e/3, 1-e]


### ACTUAL PROGRAM BELOW
def main() :
    args = InputOutput.parseArgs("AlphaPeel")
    pedigree = Pedigree.Pedigree() 
    InputOutput.readInPedigreeFromInputs(pedigree, args, haps = True)

    # Do really simple phasing.
    for ind in pedigree:
        HaplotypeOperations.setup_individual(ind)

    for ind in pedigree:
        HaplotypeOperations.fillFromParents(ind)
    # Turn really simple phasing into genotype probabilities.

    peelingInfo = PeelingInfo.createPeelingInfo(pedigree, args, noPenetrance = True)
    # simpleFill(pedigree, peelingInfo)

    # Run the segregation estimation step without anterior information (hack: anterior values all set to .25)

    estimateRecombinations(pedigree, peelingInfo, args)

    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.recomb, args.out + ".recomb", digits = 6)
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.recomb_mat, args.out + ".recomb_mat", digits = 6)
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.recomb_pat, args.out + ".recomb_pat", digits = 6)

    PeelingIO.writeGenotypes(pedigree, genoProbFunc = peelingInfo.getGenoProbs)

    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.segregation, args.out + ".seg")
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.pointSeg, args.out + ".pointSeg")

    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.penetrance, args.out + ".penetrance")
    InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.anterior, args.out + ".anterior")


    # PeelingIO.fullOutput(pedigree, peelingInfo, args)

if __name__ == "__main__":
    main()