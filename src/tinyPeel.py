import numpy as np
from numba import jit, jitclass, float32, int32, int64, optional

from tinyhouse import Pedigree
from tinyhouse import ProbMath
from tinyhouse import InputOutput 

from Peeling import Peeling
from Peeling import PeelingIO

import concurrent.futures
from itertools import repeat



def updateFamily(family, peelingInfo):
    sire = family.sire.idn
    dam = family.dam.idn
    fam = family.idn

    peelingInfo.posterior[sire,:,:] = peelingInfo.posterior[sire,:,:] / peelingInfo.posteriorSire[fam,:,:] * peelingInfo.posteriorSire_new[fam,:,:]
    peelingInfo.posterior[dam,:,:] = peelingInfo.posterior[dam,:,:] / peelingInfo.posteriorDam[fam,:,:] * peelingInfo.posteriorDam_new[fam,:,:]

    peelingInfo.posterior[sire,:,:] = peelingInfo.posterior[sire,:,:]/np.sum(peelingInfo.posterior[sire,:,:], 0)
    peelingInfo.posterior[dam,:,:] = peelingInfo.posterior[dam,:,:] / np.sum(peelingInfo.posterior[dam,:,:], 0)

    peelingInfo.posteriorSire[fam,:,:]  = peelingInfo.posteriorSire_new[fam,:,:] 
    peelingInfo.posteriorDam[fam,:,:]  = peelingInfo.posteriorDam_new[fam,:,:] 

def updateFamilyNull(family, peelingInfo):
    sire = family.sire.idn
    dam = family.dam.idn
    fam = family.idn

    peelingInfo.posterior[sire,:,:] = 1
    peelingInfo.posterior[dam,:,:] = 1

    peelingInfo.posterior[sire,:,:] = 1
    peelingInfo.posterior[dam,:,:] = 1

    peelingInfo.posteriorSire[fam,:,:]  = 1
    peelingInfo.posteriorDam[fam,:,:]  = 1

def updatePosterior(pedigree, peelingInfo, sires, dams) :

    if pedigree.mapSireToFamilies is None:
        pedigree.setupFamilyMap()

    for sire in sires:
        updateSire(sire, pedigree.mapSireToFamilies[sire], peelingInfo)

    for dam in dams:
        updateDam(dam, pedigree.mapDamToFamilies[dam], peelingInfo)


def updateSire(sire, famList, peelingInfo) :
    peelingInfo.posterior[sire,:,:] = 0
    for famId in famList:
        log_update = np.log(peelingInfo.posteriorSire_new[famId,:,:])
        peelingInfo.posterior[sire,:,:] += log_update
        peelingInfo.posteriorSire[famId,:,:] = -log_update

    for famId in famList:
        peelingInfo.posteriorSire[famId,:,:] += peelingInfo.posterior[sire,:,:]

    #Rescale values.
    peelingInfo.posterior[sire,:,:] = Peeling.expNorm1(peelingInfo.posterior[sire,:,:])
    peelingInfo.posterior[sire,:,:] /= np.sum(peelingInfo.posterior[sire,:,:], 0)

    for famId in famList:
        peelingInfo.posteriorSire[famId,:,:] = Peeling.expNorm1(peelingInfo.posteriorSire[famId,:,:])
        peelingInfo.posteriorSire[famId,:,:]  /= np.sum(peelingInfo.posteriorSire[famId,:,:], 0)
        
def updateDam(dam, famList, peelingInfo) :
    peelingInfo.posterior[dam,:,:] = 0

    for famId in famList:
        log_update = np.log(peelingInfo.posteriorDam_new[famId,:,:])
        peelingInfo.posterior[dam,:,:] += log_update

        peelingInfo.posteriorDam[famId,:,:] = -log_update

    for famId in famList:
        peelingInfo.posteriorDam[famId,:,:] += peelingInfo.posterior[dam,:,:]

    peelingInfo.posterior[dam,:,:] = Peeling.expNorm1(peelingInfo.posterior[dam,:,:])
    peelingInfo.posterior[dam,:,:] /= np.sum(peelingInfo.posterior[dam,:,:], 0)
    for famId in famList:
        peelingInfo.posteriorDam[famId,:,:] = Peeling.expNorm1(peelingInfo.posteriorDam[famId,:,:])
        peelingInfo.posteriorDam[famId,:,:]  /= np.sum(peelingInfo.posteriorDam[famId,:,:], 0)

def peelingCycle(pedigree, peelingInfo, args, singleLocusMode = False) :
    nWorkers = args.maxthreads
    for families in pedigree.genFamilies:
        print("Peeling Down")
        jit_families = [family.toJit() for family in families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers) as executor:
                 results = executor.map(Peeling.peel, jit_families, repeat(Peeling.PEEL_DOWN), repeat(peelingInfo), repeat(singleLocusMode))
        else:
            for family in jit_families:
                Peeling.peel(family, Peeling.PEEL_DOWN, peelingInfo, singleLocusMode)

    for families in reversed(pedigree.genFamilies):
        print("Peeling Up")
        jit_families = [family.toJit() for family in families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers) as executor:
                 results = executor.map(Peeling.peel, jit_families, repeat(Peeling.PEEL_UP), repeat(peelingInfo), repeat(singleLocusMode))
        else:
            # for family in jit_families:
            for family in families:
                Peeling.peel(family.toJit(), Peeling.PEEL_UP, peelingInfo, singleLocusMode)

        # for fam in families:
        #     if not args.noposterior:
        #         updateFamily(fam, peelingInfo)
        #     else:
        #         updateFamilyNull(fam, peelingInfo)
        sires = set()
        dams = set()
        for family in families:
            sires.add(family.sire.idn)
            dams.add(family.dam.idn)
        updatePosterior(pedigree, peelingInfo, sires, dams)


def getLociAndDistance(snpMap, segMap):
    nSnp = len(snpMap)
    distance = np.full(nSnp, 0, dtype = np.float32)
    loci = np.full((nSnp, 2), 0, dtype = np.int64)

    #Assume snp map and segMap are sorted.
    segIndex = 0
    for i in range(nSnp) :
        pos = snpMap[i]
        while segIndex < (len(segMap)-1) and segMap[segIndex + 1] < pos : 
            segIndex += 1
        if segIndex == 0 and segMap[segIndex] > pos :
            loci[i, :] = (segIndex, segIndex)
            distance[i] = 0
        elif segIndex == (len(segMap)-1) and segMap[segIndex] <= pos :
            loci[i, :] = (segIndex, segIndex)
            distance[i] = 0
        else:
            loci[i,:] = (segIndex, segIndex + 1)
            gap = segMap[segIndex+1] - segMap[segIndex] 
            distance[i] = 1.0 - (pos - segMap[segIndex])/gap #At distance 0, only use segIndex. At distance 1, use segIndex + 1.
            print(gap, pos, segMap[segIndex], distance[i])
    return (loci, distance)


def generateSingleLocusSegregation(peelingInfo, pedigree, args):
    segInfo = None
    print(args.segfile)
    if args.segfile is not None:
        snpMap = np.array(InputOutput.readMapFile(args.mapfile, args.startsnp, args.stopsnp)[2])
        segMap = np.array(InputOutput.readMapFile(args.segmapfile)[2])

        loci, distance = getLociAndDistance(snpMap, segMap)
        start = np.min(loci)
        stop = np.max(loci) 
        
        seg = InputOutput.readInSeg(pedigree, args.segfile, start = start, stop = stop)
        loci -= start # Re-align to seg file.
        segInfo = (seg, loci, distance)
        for i in range(len(distance)):
            segLoc0 = loci[i,0]
            segLoc1 = loci[i,1]
            peelingInfo.segregation[:,:,i] = distance[i]*seg[:,:,segLoc0] + (1-distance[i])*seg[:,:,segLoc1]

    else:
        peelingInfo.segregation[:,:,:] = .25

def runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode = False):
    #Right now maf _only_ uses the penetrance so can be estimated once.
    # if args.estmaf: Peeling.updateMaf(pedigree, peelingInfo)

    for i in range(args.ncycles):
        print("Cycle ", i)
        peelingCycle(pedigree, peelingInfo, args = args, singleLocusMode = singleLocusMode)
        peelingInfo.iteration += 1
        # if args.esttransitions: Peeling.updateSeg(peelingInfo) #Option currently disabled.
        if args.esterrors: Peeling.updatePenetrance(pedigree, peelingInfo)


### ACTUAL PROGRAM BELOW
def main() :
    args = InputOutput.parseArgs("AlphaPeel")
    pedigree = Pedigree.Pedigree() 
    InputOutput.readInPedigreeFromInputs(pedigree, args)
    pedigree.setMaf()


    for family in pedigree.getFamilies():
        ids = [ind.idx for ind in family.offspring]
        if "id-3776246612717" in ids:
            print("fam id", family.idn)
    singleLocusMode = args.runtype == "single"
    if args.runtype == "multi" and args.segfile :
        print("Running in multi-locus mode, external segfile ignored")

    # if args.runtype == "multi" :
    peelingInfo = Peeling.createPeelingInfo(pedigree, args, phaseFounder = (not args.nophasefounders))

    if singleLocusMode:
        print("Generating seg estimates")
        generateSingleLocusSegregation(peelingInfo, pedigree, args)
    runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode = singleLocusMode)

    PeelingIO.writeGenotypes(pedigree, genoProbFunc = peelingInfo.getGenoProbs)
    if not args.no_params: PeelingIO.writeOutParamaters(peelingInfo)
    if not singleLocusMode and not args.no_seg: InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.segregation, args.out + ".seg")
    # InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.penetrance, args.out + ".penetrance")
    # InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.anterior, args.out + ".anterior")
    # InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.posterior, args.out + ".posterior")

    # InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorSire, args.out + ".posteriorSire")
    # InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorDam, args.out + ".posteriorDam")

    # InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorSire_new, args.out + ".posteriorSire_new")
    # InputOutput.writeFamIndexedMatrix(pedigree, peelingInfo.posteriorDam_new, args.out + ".posteriorDam_new")


if __name__ == "__main__":
    main()