import numpy as np

from .tinyhouse import Pedigree
from .tinyhouse import InputOutput

from .Peeling import Peeling
from .Peeling import PeelingIO
from .Peeling import PeelingInfo
from .Peeling import PeelingUpdates

import concurrent.futures
from itertools import repeat
import argparse


def runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode=False):
    # Right now maf _only_ uses the penetrance so can be estimated once.
    if args.est_alt_allele_prob:
        PeelingUpdates.updateMaf(pedigree, peelingInfo)

    for i in range(args.n_cycles):
        print("Cycle ", i)
        peelingCycle(pedigree, peelingInfo, args=args, singleLocusMode=singleLocusMode)
        peelingInfo.iteration += 1

        # esttransitions is been disabled.
        # if args.esttransitions:
        #     print("Estimating the transmission rate is currently a disabled option")
        # PeelingUpdates.updateSeg(peelingInfo) #Option currently disabled.

        if args.est_geno_error_prob or args.est_seq_error_prob:
            PeelingUpdates.updatePenetrance(pedigree, peelingInfo)


def peelingCycle(pedigree, peelingInfo, args, singleLocusMode=False):
    nWorkers = args.maxthreads

    for index, generation in enumerate(pedigree.generations):
        print("Peeling Down, Generation", index)
        jit_families = [family.toJit() for family in generation.families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=nWorkers
            ) as executor:
                executor.map(
                    Peeling.peel,
                    jit_families,
                    repeat(Peeling.PEEL_DOWN),
                    repeat(peelingInfo),
                    repeat(singleLocusMode),
                )
        else:
            for family in jit_families:
                Peeling.peel(family, Peeling.PEEL_DOWN, peelingInfo, singleLocusMode)

    for index, generation in enumerate(reversed(pedigree.generations)):
        print("Peeling Up, Generation", len(pedigree.generations) - index - 1)
        jit_families = [family.toJit() for family in generation.families]

        if args.maxthreads > 1:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=nWorkers
            ) as executor:
                executor.map(
                    Peeling.peel,
                    jit_families,
                    repeat(Peeling.PEEL_UP),
                    repeat(peelingInfo),
                    repeat(singleLocusMode),
                )
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


def updatePosterior(pedigree, peelingInfo, sires, dams):
    # if pedigree.mapSireToFamilies is None or pedigree.mapDamToFamilies is None:
    #     pedigree.setupFamilyMap()

    for sire in sires:
        updateSire(sire, peelingInfo)

    for dam in dams:
        updateDam(dam, peelingInfo)


def updateSire(sire, peelingInfo):
    famList = [fam.idn for fam in sire.families]
    sire = sire.idn
    peelingInfo.posterior[sire, :, :] = 0
    for famId in famList:
        log_update = np.log(peelingInfo.posteriorSire_new[famId, :, :])
        peelingInfo.posterior[sire, :, :] += log_update
        peelingInfo.posteriorSire_minusFam[famId, :, :] = -log_update

    for famId in famList:
        peelingInfo.posteriorSire_minusFam[famId, :, :] += peelingInfo.posterior[
            sire, :, :
        ]

    # Rescale values.
    peelingInfo.posterior[sire, :, :] = Peeling.expNorm1D(
        peelingInfo.posterior[sire, :, :]
    )
    peelingInfo.posterior[sire, :, :] /= np.sum(peelingInfo.posterior[sire, :, :], 0)

    for famId in famList:
        peelingInfo.posteriorSire_minusFam[famId, :, :] = Peeling.expNorm1D(
            peelingInfo.posteriorSire_minusFam[famId, :, :]
        )
        peelingInfo.posteriorSire_minusFam[famId, :, :] /= np.sum(
            peelingInfo.posteriorSire_minusFam[famId, :, :], 0
        )


def updateDam(dam, peelingInfo):
    famList = [fam.idn for fam in dam.families]
    dam = dam.idn
    peelingInfo.posterior[dam, :, :] = 0
    for famId in famList:
        log_update = np.log(peelingInfo.posteriorDam_new[famId, :, :])
        peelingInfo.posterior[dam, :, :] += log_update
        peelingInfo.posteriorDam_minusFam[famId, :, :] = -log_update

    for famId in famList:
        peelingInfo.posteriorDam_minusFam[famId, :, :] += peelingInfo.posterior[
            dam, :, :
        ]

    peelingInfo.posterior[dam, :, :] = Peeling.expNorm1D(
        peelingInfo.posterior[dam, :, :]
    )
    peelingInfo.posterior[dam, :, :] /= np.sum(peelingInfo.posterior[dam, :, :], 0)
    for famId in famList:
        peelingInfo.posteriorDam_minusFam[famId, :, :] = Peeling.expNorm1D(
            peelingInfo.posteriorDam_minusFam[famId, :, :]
        )
        peelingInfo.posteriorDam_minusFam[famId, :, :] /= np.sum(
            peelingInfo.posteriorDam_minusFam[famId, :, :], 0
        )


def getLociAndDistance(snpMap, segMap):
    nSnp = len(snpMap)
    distance = np.full(nSnp, 0, dtype=np.float32)
    loci = np.full((nSnp, 2), 0, dtype=np.int64)

    # Assume snp map and segMap are sorted.
    segIndex = 0
    for i in range(nSnp):
        pos = snpMap[i]
        # Move along the segMap until we reach a point where we are either at the end of the map, or where the next seg marker occurs after the position in the genotype file.
        # This assumes sorting pretty heavily. Alternative would be to find the neighboring positions in the seg file for each marker in the genotype file.
        while segIndex < (len(segMap) - 1) and segMap[segIndex + 1] < pos:
            segIndex += 1

        # Now that positions are known, choose the neighboring markers and the distance to those markers.
        # First two if statements handle the begining and ends of the chromosome.
        if segIndex == 0 and segMap[segIndex] > pos:
            loci[i, :] = (segIndex, segIndex)
            distance[i] = 0
        elif segIndex == (len(segMap) - 1) and segMap[segIndex] <= pos:
            loci[i, :] = (segIndex, segIndex)
            distance[i] = 0
        else:
            loci[i, :] = (segIndex, segIndex + 1)
            gap = segMap[segIndex + 1] - segMap[segIndex]
            distance[i] = (
                1.0 - (pos - segMap[segIndex]) / gap
            )  # At distance 0, only use segIndex. At distance 1, use segIndex + 1.
    return (loci, distance)


def generateSingleLocusSegregation(peelingInfo, pedigree, args):
    if args.segfile is not None:
        # This just gets the locations in the map files.
        snpMap = np.array(
            InputOutput.readMapFile(args.map_file, args.startsnp, args.stopsnp)[2]
        )
        segMap = np.array(InputOutput.readMapFile(args.seg_map_file)[2])

        loci, distance = getLociAndDistance(snpMap, segMap)
        start = np.min(loci)
        stop = np.max(loci)

        seg = InputOutput.readInSeg(pedigree, args.seg_file, start=start, stop=stop)
        loci -= start  # Re-align to seg file.
        for i in range(len(distance)):
            segLoc0 = loci[i, 0]
            segLoc1 = loci[i, 1]
            peelingInfo.segregation[:, :, i] = (
                distance[i] * seg[:, :, segLoc0]
                + (1 - distance[i]) * seg[:, :, segLoc1]
            )

    else:
        peelingInfo.segregation[:, :, :] = 0.25


def get_probability_options():
    parse_dictionary = dict()
    parse_dictionary["geno_error_prob"] = lambda parser: parser.add_argument(
        "-geno_error_prob",
        default=0.0001,
        required=False,
        type=float,
        help="Genotyping error rate. Default: 0.0001.",
    )
    parse_dictionary["seq_error_prob"] = lambda parser: parser.add_argument(
        "-seq_error_prob",
        default=0.001,
        required=False,
        type=float,
        help="Sequencing error rate. Default: 0.001.",
    )
    parse_dictionary["recombination"] = lambda parser: parser.add_argument(
        "-recomb",
        default=1,
        required=False,
        type=float,
        help="Recombination rate per chromosome (in Morgan). Default: 1.",
    )
    return parse_dictionary


def get_input_options():
    parse_dictionary = dict()
    parse_dictionary["bfile"] = lambda parser: parser.add_argument(
        "-plink_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A file in plink (binary) format. Only stable on Linux.",
    )
    parse_dictionary["genotypes"] = lambda parser: parser.add_argument(
        "-geno_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A file in AlphaGenes format.",
    )
    parse_dictionary["reference"] = lambda parser: parser.add_argument(
        "-reference",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A haplotype reference panel in AlphaGenes format.",
    )
    parse_dictionary["seqfile"] = lambda parser: parser.add_argument(
        "-seq_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A sequence data file.",
    )
    parse_dictionary["pedigree"] = lambda parser: parser.add_argument(
        "-ped_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A pedigree file in AlphaGenes format.",
    )
    parse_dictionary["phasefile"] = lambda parser: parser.add_argument(
        "-hap_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="A haplotype file in AlphaGenes format.",
    )
    parse_dictionary["startsnp"] = lambda parser: parser.add_argument(
        "-start_snp",
        default=None,
        required=False,
        type=int,
        help="The first marker to consider. The first marker in the file is marker '1'. Default: 1.",
    )
    parse_dictionary["stopsnp"] = lambda parser: parser.add_argument(
        "-stop_snp",
        default=None,
        required=False,
        type=int,
        help="The last marker to consider. Default: all markers considered.",
    )
    parse_dictionary["seed"] = lambda parser: parser.add_argument(
        "-seed",
        default=None,
        required=False,
        type=int,
        help="A random seed to use for debugging.",
    )

    return parse_dictionary


def get_output_options():
    parse_dictionary = dict()

    parse_dictionary["writekey"] = lambda parser: parser.add_argument(
        "-out_id_order",
        default="id",
        required=False,
        type=str,
        help='Determines the order in which individuals are ordered in the output file based on their order in the corresponding input file. Animals not in the input file are placed at the end of the file and sorted in alphanumeric order. These animals can be surpressed with the "-onlykeyed" option. Options: id, pedigree, genotypes, sequence, segregation. Defualt: id.',
    )
    parse_dictionary["onlykeyed"] = lambda parser: parser.add_argument(
        "-out_id_only",
        action="store_true",
        required=False,
        help='Flag to surpress the animals who are not present in the file used with -outputkey. Also surpresses "dummy" animals.',
    )
    parse_dictionary["iothreads"] = lambda parser: parser.add_argument(
        "-n_io_threads",
        default=1,
        required=False,
        type=int,
        help="Number of threads to use for io. Default: 1.",
    )

    return parse_dictionary


def get_multithread_options():
    parse_dictionary = dict()
    parse_dictionary["iothreads"] = lambda parser: parser.add_argument(
        "-n_io_threads",
        default=1,
        required=False,
        type=int,
        help="Number of threads to use for input and output. Default: 1.",
    )
    parse_dictionary["maxthreads"] = lambda parser: parser.add_argument(
        "-n_threads",
        default=1,
        required=False,
        type=int,
        help="Maximum number of threads to use for analysis. Default: 1.",
    )
    return parse_dictionary


# ACTUAL PROGRAM BELOW


def getArgs():
    parser = argparse.ArgumentParser(description="")
    core_parser = parser.add_argument_group("Core arguments")
    core_parser.add_argument(
        "-out_file", required=True, type=str, help="The output file prefix."
    )

    core_peeling_parser = parser.add_argument_group("Mandatory peeling arguments")
    core_peeling_parser.add_argument(
        "-method",
        default=None,
        required=False,
        type=str,
        help='Program run type. Either "single" or "multi".',
    )

    # Input options
    input_parser = parser.add_argument_group("Input Options")
    InputOutput.add_arguments_from_dictionary(
        input_parser,
        get_input_options(),
        options=[
            "bfile",
            "genotypes",
            "phasefile",
            "seqfile",
            "pedigree",
            "startsnp",
            "stopsnp",
        ],
    )

    # Output options
    output_parser = parser.add_argument_group("Output Options")

    output_parser.add_argument(
        "-no_dosage",
        action="store_true",
        required=False,
        help="Flag to suppress the output of the dosage file.",
    )
    output_parser.add_argument(
        "-seg_prob",
        action="store_true",
        required=False,
        help="Flag to enable writing out the segregation probabilities.",
    )
    output_parser.add_argument(
        "-no_params",
        action="store_true",
        required=False,
        help="Flag to suppress writing the model parameter files.",
    )
    output_parser.add_argument(
        "-geno_prob",
        action="store_true",
        required=False,
        help="Flag to enable writing out the genotype probabilities.",
    )
    output_parser.add_argument(
        "-phased_geno_prob",
        action="store_true",
        required=False,
        help="Flag to enable writing out the phased genotype probabilities.",
    )
    output_parser.add_argument(
        "-geno_threshold",
        default=None,
        required=False,
        type=float,
        nargs="*",
        help="Genotype calling threshold(s). Multiple space separated values allowed.\
        Value less than 1/3 will be replaced by 1/3.",
    )
    output_parser.add_argument(
        "-hap_threshold",
        default=None,
        required=False,
        type=float,
        nargs="*",
        help="Haplotype calling threshold(s). Multiple space separated values allowed.\
        Value less than 1/2 will be replaced by 1/2.",
    )
    output_parser.add_argument(
        "-geno",
        action="store_true",
        required=False,
        help="Flag to call and write out the genotypes.",
    )
    output_parser.add_argument(
        "-binary_call_files",
        action="store_true",
        required=False,
        help="Flag to write out the called genotype files as a binary plink output [Not yet implemented].",
    )
    output_parser.add_argument(
        "-hap",
        action="store_true",
        required=False,
        help="Flag to call and write out the haplotypes.",
    )

    InputOutput.add_arguments_from_dictionary(
        output_parser,
        get_output_options(),
        options=["writekey", "onlykeyed"],
    )

    # Multithreading
    multithread_parser = parser.add_argument_group("Multithreading Options")
    InputOutput.add_arguments_from_dictionary(
        multithread_parser,
        get_multithread_options(),
        options=["iothreads", "maxthreads"],
    )

    peeling_parser = parser.add_argument_group("Optional peeling arguments")
    peeling_parser.add_argument(
        "-n_cycles",
        default=5,
        required=False,
        type=int,
        help="Number of peeling cycles. Default: 5.",
    )
    peeling_parser.add_argument(
        "-rec_length",
        default=1.0,
        required=False,
        type=float,
        help="Estimated recombination length of the chromosome in Morgans. [Default 1.00]",
    )
    peeling_parser.add_argument(
        "-penetrance",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help=argparse.SUPPRESS,
    )  # help='An optional external penetrance file. This will overwrite the default penetrance values.')
    InputOutput.add_arguments_from_dictionary(
        peeling_parser,
        get_probability_options(),
        options=["geno_error_prob", "seq_error_prob"],
    )

    peeling_control_parser = parser.add_argument_group("Peeling control arguments")
    peeling_control_parser.add_argument(
        "-est_geno_error_prob",
        action="store_true",
        required=False,
        help="Flag to re-estimate the genotyping error rates after each peeling cycle.",
    )
    peeling_control_parser.add_argument(
        "-est_seq_error_prob",
        action="store_true",
        required=False,
        help="Flag to re-estimate the sequencing error rates after each peeling cycle.",
    )
    peeling_control_parser.add_argument(
        "-est_alt_allele_prob",
        action="store_true",
        required=False,
        help="Flag to re-estimate the alternative allele probabilities after each peeling cycle.",
    )
    peeling_control_parser.add_argument(
        "-nophasefounders",
        action="store_true",
        required=False,
        help="A flag phase a heterozygous allele in one of the founders (if such an allele can be found).",
    )
    peeling_control_parser.add_argument(
        "-sex_chrom",
        action="store_true",
        required=False,
        help="A flag to indicate that input data is for a sex chromosome. Sex needs to be given in the pedigree file. This is currently an experimental option.",
    )

    singleLocus_parser = parser.add_argument_group("Hybrid peeling arguments")
    singleLocus_parser.add_argument(
        "-map_file",
        default=None,
        required=False,
        type=str,
        help="a map file for genotype data.",
    )
    singleLocus_parser.add_argument(
        "-seg_map_file",
        default=None,
        required=False,
        type=str,
        help="a map file for the segregation estimates for hybrid peeling.",
    )
    singleLocus_parser.add_argument(
        "-seg_file",
        default=None,
        required=False,
        type=str,
        help="A segregation file for hybrid peeling.",
    )
    # singleLocus_parser.add_argument('-blocksize',default=100, required=False, type=int, help='The number of markers to impute at once. This changes the memory requirements of the program.')

    return InputOutput.parseArgs("AlphaPeel", parser)


def main():
    args = getArgs()
    pedigree = Pedigree.Pedigree()
    InputOutput.readInPedigreeFromInputs(pedigree, args)

    singleLocusMode = args.method == "single"
    if args.method == "multi" and args.segfile:
        print("Running in multi-locus mode, external segfile ignored")

    peelingInfo = PeelingInfo.createPeelingInfo(
        pedigree, args, phaseFounder=(not args.nophasefounders)
    )

    if singleLocusMode:
        print("Generating seg estimates")
        generateSingleLocusSegregation(peelingInfo, pedigree, args)
    runPeelingCycles(pedigree, peelingInfo, args, singleLocusMode=singleLocusMode)

    PeelingIO.writeGenotypes(
        pedigree,
        genoProbFunc=peelingInfo.getGenoProbs,
        isSexChrom=peelingInfo.isSexChrom,
    )
    if not args.no_params:
        PeelingIO.writeOutParamaters(peelingInfo)
    if not singleLocusMode and args.seg_prob:
        InputOutput.writeIdnIndexedMatrix(
            pedigree, peelingInfo.segregation, args.out_file + ".seg_prob.txt"
        )


if __name__ == "__main__":
    main()
