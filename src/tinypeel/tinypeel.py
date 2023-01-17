import numpy as np
from numba import jit, float32, int32, int64, optional
from numba.experimental import jitclass

from .tinyhouse import Pedigree
from .tinyhouse import ProbMath
from .tinyhouse import InputOutput 

from .Peeling import Peeling
from .Peeling import PeelingIO
from .Peeling import PeelingInfo
from .Peeling import PeelingUpdates

import concurrent.futures
from itertools import repeat
import argparse

def runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode = False):
    #Right now maf _only_ uses the penetrance so can be estimated once.
    if args.estmaf: PeelingUpdates.updateMaf(pedigree, peelingInfo)

    for i in range(args.ncycles):
        print("Cycle ", i)
        peelingCycle(pedigree, peelingInfo, args = args, singleLocusMode = singleLocusMode)
        peelingInfo.iteration += 1
        
        # esttransitions is been disabled.
        # if args.esttransitions: 
        #     print("Estimating the transmission rate is currently a disabled option")
            # PeelingUpdates.updateSeg(peelingInfo) #Option currently disabled.
            
        if args.esterrors: 
            PeelingUpdates.updatePenetrance(pedigree, peelingInfo)

def peelingCycle(pedigree, peelingInfo, args, singleLocusMode = False) :
    nWorkers = args.maxthreads
    
    for index, generation in enumerate(pedigree.generations):
        print("Peeling Down, Generation", index)
        jit_families = [family.toJit() for family in generation.families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers) as executor:
                 results = executor.map(Peeling.peel, jit_families, repeat(Peeling.PEEL_DOWN), repeat(peelingInfo), repeat(singleLocusMode))
        else:
            for family in jit_families:
                Peeling.peel(family, Peeling.PEEL_DOWN, peelingInfo, singleLocusMode)

    for index, generation in enumerate(reversed(pedigree.generations)):
        print("Peeling Up, Generation", len(pedigree.generations) - index - 1)
        jit_families = [family.toJit() for family in generation.families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nWorkers) as executor:
                 results = executor.map(Peeling.peel, jit_families, repeat(Peeling.PEEL_UP), repeat(peelingInfo), repeat(singleLocusMode))
        else:
            for family in jit_families:
                Peeling.peel(family, Peeling.PEEL_UP, peelingInfo, singleLocusMode)

        sires = set()
        dams = set()
        for family in generation.families:
            sires.add(family.sire)
            dams.add(family.dam)
        updatePosterior(pedigree, peelingInfo, sires, dams)

# updatePosterior updates the posterior term for a specific set of sires and dams.
# The updateSire and updateDam functions perform the updates for a specific sire 
# and specific dam by including all of the information from all of their families.
# This update is currently not multithreaded. It is also currently not ideal -- 
# right now the posterior term is updated for all of the sires/dams no matter 
# whether or not they have been changed since the last update.

def updatePosterior(pedigree, peelingInfo, sires, dams) :

    # if pedigree.mapSireToFamilies is None or pedigree.mapDamToFamilies is None:
    #     pedigree.setupFamilyMap()

    for sire in sires:
        updateSire(sire, peelingInfo)

    for dam in dams:
        updateDam(dam, peelingInfo)


def updateSire(sire, peelingInfo) :

    famList = [fam.idn for fam in sire.families]
    sire = sire.idn
    peelingInfo.posterior[sire,:,:] = 0
    for famId in famList:
        log_update = np.log(peelingInfo.posteriorSire_new[famId,:,:])
        peelingInfo.posterior[sire,:,:] += log_update
        peelingInfo.posteriorSire_minusFam[famId,:,:] = -log_update

    for famId in famList:
        peelingInfo.posteriorSire_minusFam[famId,:,:] += peelingInfo.posterior[sire,:,:]

    #Rescale values.
    peelingInfo.posterior[sire,:,:] = Peeling.expNorm1D(peelingInfo.posterior[sire,:,:])
    peelingInfo.posterior[sire,:,:] /= np.sum(peelingInfo.posterior[sire,:,:], 0)

    for famId in famList:
        peelingInfo.posteriorSire_minusFam[famId,:,:] = Peeling.expNorm1D(peelingInfo.posteriorSire_minusFam[famId,:,:])
        peelingInfo.posteriorSire_minusFam[famId,:,:]  /= np.sum(peelingInfo.posteriorSire_minusFam[famId,:,:], 0)
        
def updateDam(dam, peelingInfo) :

    famList = [fam.idn for fam in dam.families]
    dam = dam.idn
    peelingInfo.posterior[dam,:,:] = 0
    for famId in famList:
        log_update = np.log(peelingInfo.posteriorDam_new[famId,:,:])
        peelingInfo.posterior[dam,:,:] += log_update
        peelingInfo.posteriorDam_minusFam[famId,:,:] = -log_update

    for famId in famList:
        peelingInfo.posteriorDam_minusFam[famId,:,:] += peelingInfo.posterior[dam,:,:]

    peelingInfo.posterior[dam,:,:] = Peeling.expNorm1D(peelingInfo.posterior[dam,:,:])
    peelingInfo.posterior[dam,:,:] /= np.sum(peelingInfo.posterior[dam,:,:], 0)
    for famId in famList:
        peelingInfo.posteriorDam_minusFam[famId,:,:] = Peeling.expNorm1D(peelingInfo.posteriorDam_minusFam[famId,:,:])
        peelingInfo.posteriorDam_minusFam[famId,:,:]  /= np.sum(peelingInfo.posteriorDam_minusFam[famId,:,:], 0)

def getLociAndDistance(snpMap, segMap):
    nSnp = len(snpMap)
    distance = np.full(nSnp, 0, dtype = np.float32)
    loci = np.full((nSnp, 2), 0, dtype = np.int64)

    # Assume snp map and segMap are sorted.
    segIndex = 0
    for i in range(nSnp) :
        pos = snpMap[i]
        # Move along the segMap until we reach a point where we are either at the end of the map, or where the next seg marker occurs after the position in the genotype file.
        # This assumes sorting pretty heavily. Alternative would be to find the neighboring positions in the seg file for each marker in the genotype file.
        while segIndex < (len(segMap)-1) and segMap[segIndex + 1] < pos : 
            segIndex += 1
        
        # Now that positions are known, choose the neighboring markers and the distance to those markers.
        # First two if statements handle the begining and ends of the chromosome.
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
    return (loci, distance)


def generateSingleLocusSegregation(peelingInfo, pedigree, args):
    segInfo = None
    if args.segfile is not None:
        # This just gets the locations in the map files.
        snpMap = np.array(InputOutput.readMapFile(args.mapfile, args.startsnp, args.stopsnp)[2])
        segMap = np.array(InputOutput.readMapFile(args.segmapfile)[2])

        loci, distance = getLociAndDistance(snpMap, segMap)
        start = np.min(loci)
        stop = np.max(loci) 
        
        seg = InputOutput.readInSeg(pedigree, args.segfile, start = start, stop = stop)
        loci -= start # Re-align to seg file.
        for i in range(len(distance)):
            segLoc0 = loci[i,0]
            segLoc1 = loci[i,1]
            peelingInfo.segregation[:,:,i] = distance[i]*seg[:,:,segLoc0] + (1-distance[i])*seg[:,:,segLoc1]

    else:
        peelingInfo.segregation[:,:,:] = .25


### ACTUAL PROGRAM BELOW


def getArgs() :
    parser = argparse.ArgumentParser(description='')
    core_parser = parser.add_argument_group("Core arguments")
    core_parser.add_argument('-out', required=True, type=str, help='The output file prefix.')
   
    core_peeling_parser = parser.add_argument_group("Mandatory peeling arguments")
    core_peeling_parser.add_argument('-runtype', default=None, required=False, type=str, help='Program run type. Either "single" or "multi".')

    # Input options
    input_parser = parser.add_argument_group("Input Options")
    InputOutput.add_arguments_from_dictionary(input_parser, InputOutput.get_input_options(), options = ["bfile", "genotypes", "phasefile", "seqfile", "pedigree", "startsnp", "stopsnp"]) 

    # Output options
    output_parser = parser.add_argument_group("Output Options")

    output_parser.add_argument('-no_dosages', action='store_true', required=False, help='Flag to suppress the dosage files.')
    output_parser.add_argument('-no_seg', action='store_true', required=False, help='Flag to suppress the segregation files (e.g. when running for chip imputation and not hybrid peeling).')
    output_parser.add_argument('-no_params', action='store_true', required=False, help='Flag to suppress writing the parameter files.')

    output_parser.add_argument('-haps', action='store_true', required=False, help='Flag to enable writing out the genotype probabilities.')
    output_parser.add_argument('-calling_threshold', default=None, required=False, type=float, nargs="*", help='Genotype calling threshold(s). Multiple space separated values allowed. Use. .3 for best guess genotype.')
    output_parser.add_argument('-binary_call_files', action='store_true', required=False, help='Flag to write out the called genotype files as a binary plink output [Not yet implemented].')
    output_parser.add_argument('-call_phase', action='store_true', required=False, help='Flag to call the phase as well as the genotypes.')
    
    InputOutput.add_arguments_from_dictionary(output_parser, InputOutput.get_output_options(), options = ["writekey", "onlykeyed"]) 


    # Multithreading
    multithread_parser = parser.add_argument_group("Multithreading Options")
    InputOutput.add_arguments_from_dictionary(multithread_parser, InputOutput.get_multithread_options(), options = ["iothreads", "maxthreads"]) 


    peeling_parser = parser.add_argument_group("Optional peeling arguments")
    peeling_parser.add_argument('-ncycles',default=5, required=False, type=int, help='Number of peeling cycles. Default: 5.')
    peeling_parser.add_argument('-length', default=1.0, required=False, type=float, help='Estimated length of the chromosome in Morgans. [Default 1.00]')
    peeling_parser.add_argument('-penetrance',   default=None, required=False, type=str, nargs="*", help=argparse.SUPPRESS) #help='An optional external penetrance file. This will overwrite the default penetrance values.')
    InputOutput.add_arguments_from_dictionary(peeling_parser, InputOutput.get_probability_options(), options = ["error", "seqerror"]) 


    peeling_control_parser = parser.add_argument_group("Peeling control arguments")
    peeling_control_parser.add_argument('-esterrors', action='store_true', required=False, help='Flag to re-estimate the genotyping error rates after each peeling cycle.')
    peeling_control_parser.add_argument('-estmaf', action='store_true', required=False, help='Flag to re-estimate the minor allele frequency after each peeling cycle.')
    peeling_control_parser.add_argument('-nophasefounders', action='store_true', required=False, help='A flag phase a heterozygous allele in one of the founders (if such an allele can be found).')
    peeling_control_parser.add_argument('-sexchrom', action='store_true', required=False, help='A flag to that this is a sex chromosome. Sex needs to be given in the pedigree file. This is currently an experimental option.')


    singleLocus_parser = parser.add_argument_group("Hybrid peeling arguments")
    singleLocus_parser.add_argument('-mapfile',default=None, required=False, type=str, help='a map file for genotype data.')
    singleLocus_parser.add_argument('-segmapfile',default=None, required=False, type=str, help='a map file for the segregation estimates for hybrid peeling.')
    singleLocus_parser.add_argument('-segfile',default=None, required=False, type=str, help='A segregation file for hybrid peeling.')
    # singleLocus_parser.add_argument('-blocksize',default=100, required=False, type=int, help='The number of markers to impute at once. This changes the memory requirements of the program.')




    return InputOutput.parseArgs("AlphaPeel", parser)





def main() :
    args = getArgs()
    pedigree = Pedigree.Pedigree() 
    InputOutput.readInPedigreeFromInputs(pedigree, args)

    singleLocusMode = args.runtype == "single"
    if args.runtype == "multi" and args.segfile :
        print("Running in multi-locus mode, external segfile ignored")

    peelingInfo = PeelingInfo.createPeelingInfo(pedigree, args, phaseFounder = (not args.nophasefounders))

    if singleLocusMode:
        print("Generating seg estimates")
        generateSingleLocusSegregation(peelingInfo, pedigree, args)
    runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode = singleLocusMode)

    PeelingIO.writeGenotypes(pedigree, genoProbFunc = peelingInfo.getGenoProbs)
    if not args.no_params: PeelingIO.writeOutParamaters(peelingInfo)
    if not singleLocusMode and not args.no_seg: InputOutput.writeIdnIndexedMatrix(pedigree, peelingInfo.segregation, args.out + ".seg")


if __name__ == "__main__":
    main()
