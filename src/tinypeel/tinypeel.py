"""Entry point for the AlphaPeel program."""

import sys
import warnings
import concurrent.futures
import argparse
from itertools import repeat

import numpy as np

from .tinyhouse.Pedigree import Pedigree
from .tinyhouse.InputOutput import (
    print_boilerplate,
    readMapFile,
    readInSeg,
    add_arguments_from_dictionary,
    parseArgs,
    readInPedigreeFromInputs,
)

from .peeling.peeling import (
    peel_down,
    peel_up,
    update_posterior,
)
from .peeling.peeling_io import write_requested_outputs
from .peeling.peeling_info_module import (
    create_peeling_info,
    create_locus_block_peeling_infos,
    PeelingCycleContext,
)
from .peeling.peeling_updates import (
    prepare_alternative_allele_probabilities,
    update_penetrance,
    update_maf_after_peeling,
    update_pheno_penetrance,
)
from .peeling.version import version

version_version = version

ALPHAPEEL_ARGUMENT_ALIASES = {
    "pedigree": "ped_file",
    "genotypes": "geno_file",
    "seqfile": "seq_file",
    "phenotype": "pheno_file",
    "phenoPenetrance": "pheno_penetrance_prob_file",
    "phasefile": "hap_file",
    "writekey": "out_id_order",
    "onlykeyed": "out_id_only",
    "maxthreads": "n_thread_fam",
    "segfile": "seg_file",
}


def run_peeling_cycles(pedigree, peeling_info, args, single_locus_mode=False):
    """Set up and run the configured peeling cycles.

    The setup prepares alternative allele probabilities for each metafounder.
    These priors either come from ``-alt_allele_prob_file`` or default to 0.5,
    and ``-est_start_alt_allele_prob`` can replace them with penetrance-based
    estimates before cycling. Each cycle then peels down and up through the
    pedigree generations and updates requested model parameters.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :param single_locus_mode: whether method is single locus or not, defaults to False
    :type single_locus_mode: bool, optional
    :return: None. The function modifies the peeling_info and pedigree object in place
    """
    # Initial MAF estimates depend only on penetrance, so they can be prepared once.
    prepare_alternative_allele_probabilities(pedigree, peeling_info, args)

    locus_thread_blocks = None
    if args.n_thread_loci > 1:
        locus_thread_blocks = create_locus_block_peeling_infos(
            peeling_info, args.n_thread_loci
        )
        if len(locus_thread_blocks) <= 1:
            locus_thread_blocks = None

    jit_generations = None
    if args.n_cycle > 0:
        jit_generations = get_jit_families_by_generation(pedigree)

    cycle_context = PeelingCycleContext(
        pedigree=pedigree,
        peeling_info=peeling_info,
        n_fam_threads=args.maxthreads,
        single_locus_mode=single_locus_mode,
        jit_generations=jit_generations,
        locus_thread_blocks=locus_thread_blocks,
    )

    for i in range(args.n_cycle):
        print("Cycle ", i)
        peeling_cycle(cycle_context)
        peeling_info.iteration += 1
        update_estimated_parameters(pedigree, peeling_info, args)


def update_estimated_parameters(pedigree, peeling_info, args):
    """Update requested estimated parameters after one peeling cycle."""

    if args.est_geno_error_prob or args.est_seq_error_prob:
        update_penetrance(pedigree, peeling_info, args)
    if args.est_pheno_penetrance_prob:
        update_estimated_phenotype_penetrance(pedigree, peeling_info, args)
    if args.est_alt_allele_prob:
        print("Updating Alternative Allele Frequencies")
        update_maf_after_peeling(pedigree, peeling_info)


def update_estimated_phenotype_penetrance(pedigree, peeling_info, args):
    """Update phenotype penetrance when required inputs are available."""

    if args.phenoPenetrance is None or args.phenotype is None:
        warnings.warn(
            "Both -pheno_penetrance_prob_file and -pheno_file are "
            "required to update the phenotype penetrance probabilities. "
            "Skipping update."
        )
        return

    print("Updating Phenotype Penetrance")
    update_pheno_penetrance(pedigree, peeling_info)


def get_jit_families_by_generation(pedigree):
    """Build reusable jit family containers for each generation."""

    return [
        [family.toJit() for family in generation.families]
        for generation in pedigree.generations
    ]


def peeling_cycle(context):
    """Run one peeling cycle, first down the pedigree and then back up.

    :param context: prepared state for this peeling cycle
    :type context: PeelingCycleContext
    :return: None. The function modifies the peeling_info and pedigree object in place
    """
    pedigree = context.pedigree

    if context.jit_generations is None:
        context.jit_generations = get_jit_families_by_generation(pedigree)

    for index, jit_families in enumerate(context.jit_generations):
        peel_down_generation(context, index, jit_families)

    for index, generation in enumerate(reversed(pedigree.generations)):
        generation_index = pedigree.nGenerations - index - 1
        peel_up_generation(context, generation_index, generation)


def peel_down_generation(context, generation_index, jit_families):
    """Run peel-down for one generation."""

    print("Peeling Down, Generation", generation_index)

    if context.n_fam_threads > 1:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=context.n_fam_threads
        ) as executor:
            executor.map(
                peel_down,
                jit_families,
                repeat(context.peeling_info),
                repeat(context.single_locus_mode),
                repeat(context.locus_thread_blocks),
            )
        return

    for family in jit_families:
        peel_down(
            family,
            context.peeling_info,
            context.single_locus_mode,
            context.locus_thread_blocks,
        )


def peel_up_generation(context, generation_index, generation):
    """Run peel-up for one generation and update affected parent posteriors."""

    print("Peeling Up, Generation", generation_index)
    jit_families = context.jit_generations[generation_index]

    if context.n_fam_threads > 1:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=context.n_fam_threads
        ) as executor:
            executor.map(
                peel_up,
                jit_families,
                repeat(context.peeling_info),
                repeat(context.locus_thread_blocks),
            )
    else:
        for family in jit_families:
            peel_up(family, context.peeling_info, context.locus_thread_blocks)

    sires = set()
    dams = set()
    for family in generation.families:
        sires.add(family.sire)
        dams.add(family.dam)
    update_posterior(context.peeling_info, sires, dams)


def get_loci_and_distance(snp_map, seg_map):
    """Return the bracketing segregation loci and interpolation distance for each SNP.

    ``snp_map`` contains the genotype/SNP positions and ``seg_map`` contains
    the segregation-marker positions. Both must already be sorted. The returned
    ``loci`` values identify the two neighbouring segregation markers for each
    SNP, while ``distance`` is the interpolation weight between those markers.

    :param snp_map: matrix for the position of the SNPs across and on the chromosome.
    :type snp_map: 1D numpy array with length of n_loci
    :param seg_map: matrix for the position of the segregation markers across and on the chromosome.
    :type seg_map: 1D numpy array with length of n_loci
    :return: loci and distance for each SNP
    :rtype: tuple of (loci, distance)
    """
    n_snp = len(snp_map)
    distance = np.full(n_snp, 0, dtype=np.float32)
    loci = np.full((n_snp, 2), 0, dtype=np.uint32)

    # Both maps must be sorted so the segregation index can advance monotonically.
    seg_index = 0
    for i in range(n_snp):
        pos = snp_map[i]
        # Find the closest marker at or before the SNP position.
        while seg_index < (len(seg_map) - 1) and seg_map[seg_index + 1] < pos:
            seg_index += 1

        # Clamp positions outside the segregation map to the nearest boundary.
        if seg_index == 0 and seg_map[seg_index] > pos:
            loci[i, 0] = seg_index
            loci[i, 1] = seg_index
            distance[i] = 0
        elif seg_index == (len(seg_map) - 1) and seg_map[seg_index] <= pos:
            loci[i, 0] = seg_index
            loci[i, 1] = seg_index
            distance[i] = 0
        else:
            loci[i, 0] = seg_index
            loci[i, 1] = seg_index + 1
            gap = seg_map[seg_index + 1] - seg_map[seg_index]
            distance[i] = 1.0 - (pos - seg_map[seg_index]) / gap
    return (loci, distance)


def generate_single_locus_segregation(peeling_info, pedigree, args):
    """Populate single-locus segregation probabilities from a segregation map.

    When ``-seg_file`` is provided, ``peeling_info.positions`` supplies SNP
    positions from the genotype map and ``-seg_map_file`` supplies segregation
    marker positions. The function reads only the required ``start:stop`` window
    from the segregation file, shifts the matched loci back to that local window,
    and interpolates each SNP from its two neighbouring segregation markers.
    Otherwise the default uniform probabilities remain in place.

    :param peeling_info: Peeling information container
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :return: None. The function modifies the peeling_info object in place
    """
    if args.segfile is not None:
        seg, loci, distance = read_single_locus_segregation_inputs(
            peeling_info, pedigree, args
        )
        interpolate_single_locus_segregation(
            peeling_info.segregation, seg, loci, distance
        )


def read_single_locus_segregation_inputs(peeling_info, pedigree, args):
    """Read segregation inputs and map SNP loci into the local segregation window."""

    seg_map = np.array(readMapFile(args.seg_map_file)[2])
    loci, distance = get_loci_and_distance(peeling_info.positions, seg_map)
    start = np.min(loci)
    stop = np.max(loci)
    seg = readInSeg(pedigree, args.seg_file, start=start, stop=stop)
    # Re-align absolute segregation-map indices to the window read from seg_file.
    loci -= start
    return seg, loci, distance


def interpolate_single_locus_segregation(segregation, seg, loci, distance):
    """Interpolate each SNP's segregation probabilities from bracketing markers."""

    for i, value in enumerate(distance):
        seg_loc0 = loci[i, 0]
        seg_loc1 = loci[i, 1]
        segregation_at_locus = segregation[:, :, i]
        seg0 = seg[:, :, seg_loc0]
        seg1 = seg[:, :, seg_loc1]
        segregation_at_locus[:, :] = value * seg0 + (1 - value) * seg1


def get_probability_options():
    """Collects potential user inputs for genotype error rate and sequencing error rate,
    otherwise the default is 0.0001 and 0.001 respectively.

    :return: A dictionary with the options for the genotype and sequencing error rates.
    :rtype: dict
    """
    parse_dictionary = {}
    parse_dictionary["mut_prob"] = lambda parser: parser.add_argument(
        "-mut_prob",
        default=1e-8,
        required=False,
        type=float,
        help="Mutation probability. Default: 1e-8.",
    )
    parse_dictionary["geno_error_prob"] = lambda parser: parser.add_argument(
        "-geno_error_prob",
        default=0.0001,
        required=False,
        type=float,
        help="Genotype error probability. Default: 0.0001.",
    )
    parse_dictionary["seq_error_prob"] = lambda parser: parser.add_argument(
        "-seq_error_prob",
        default=0.001,
        required=False,
        type=float,
        help="Sequencing error probability. Must not be 0. Default: 0.001.",
    )
    return parse_dictionary


def get_input_options():
    """Collects the input options of the program as a dictionary. The options are:
    -geno_file: genotypes
    -pheno_file: phenotype
    -seq_file: seqfile
    -ped_file: pedigree
    -hap_file: phasefile
    -alt_allele_prob_file: alt_allele_prob_file
    -pheno_penetrance_prob_file: pheno_penetrance_prob_file
    -start_snp: startsnp
    -stop_snp: stopsnp
    -main_metafounder: main_metafounder
    -seed: seed (for debugging)

    :return: the options for the input files and parameters.
    :rtype: dict
    """
    parse_dictionary = {}
    parse_dictionary["pedigree"] = lambda parser: parser.add_argument(
        "-ped_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Pedigree file(s) in (see format details in the docs).",
    )
    parse_dictionary["genotypes"] = lambda parser: parser.add_argument(
        "-geno_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Genotype file(s) (see format details in the docs).",
    )
    parse_dictionary["seqfile"] = lambda parser: parser.add_argument(
        "-seq_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Sequence allele read count file(s) (see format details in the docs).",
    )

    parse_dictionary["startsnp"] = lambda parser: parser.add_argument(
        "-start_snp",
        default=None,
        required=False,
        type=int,
        help="The first locus to consider. Counting starts at 1. Default: 1.",
    )
    parse_dictionary["stopsnp"] = lambda parser: parser.add_argument(
        "-stop_snp",
        default=None,
        required=False,
        type=int,
        help="The last locus to consider. Default: all loci in input files",
    )
    parse_dictionary["alt_allele_prob_file"] = lambda parser: parser.add_argument(
        "-alt_allele_prob_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Alternative allele probability file (see format details in the docs). "
        "Default: 0.5 for each locus.",
    )
    parse_dictionary["main_metafounder"] = lambda parser: parser.add_argument(
        "-main_metafounder",
        default="MF_1",
        required=False,
        type=str,
        help="ID used to represent the base population of the pedigree "
        "(metafounder / unknown parent group) (see format details in the docs). "
        "Default: MF_1.",
    )
    parse_dictionary["phenotype"] = lambda parser: parser.add_argument(
        "-pheno_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Phenotype file(s) (see format details in the docs).",
    )
    parse_dictionary["pheno_penetrance_prob_file"] = lambda parser: parser.add_argument(
        "-pheno_penetrance_prob_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Phenotype penetrance probability file (see format details in the docs).",
    )
    parse_dictionary["phasefile"] = lambda parser: parser.add_argument(
        "-hap_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Haplotype file (see format details in the docs).",
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
    """Collects the optional output options of the program as a dictionary. The options are:
    -out_id_order: writekey (write in the same order as pedigree input)
    -out_id_only: onlykeyed (only include individuals from the pedigree input)

    :return: the options for the output files and parameters.
    :rtype: dict
    """
    parse_dictionary = {}

    parse_dictionary["writekey"] = lambda parser: parser.add_argument(
        "-out_id_order",
        default="id",
        required=False,
        type=str,
        help="Determines the order in which individuals are ordered in the "
        "output file based on their order in the corresponding input file. "
        "Individuals not in the input file are placed at the end of the file and "
        "sorted in alphanumeric order. "
        'These inividuals can be surpressed with the "-out_id_only" option. '
        "Options: id, pedigree, genotypes, sequence, segregation. Defualt: id.",
    )
    parse_dictionary["onlykeyed"] = lambda parser: parser.add_argument(
        "-out_id_only",
        action="store_true",
        required=False,
        help="Suppress output for individuals not present in the file "
        "specified with -out_id_order. "
        'It also suppresses "dummy" individuals.',
    )
    parse_dictionary["out_digits"] = lambda parser: parser.add_argument(
        "-out_digits",
        default=4,
        type=int,
        required=False,
        help="Specify the number of digits to round the outputs. "
        "Does not apply to outputs from ``alt_allele_prob``, ``geno_error_prob``, "
        "``seq_error_prob``, and ``pheno_penetrance``. Default: 4.",
    )

    return parse_dictionary


def get_multithread_options():
    """Collects the optional multithread options of the program as a dictionary. The option is:
    -n_thread_fam: maxthreads

    :return: the option for the multithreading parameters.
    :rtype: dict
    """
    parse_dictionary = {}
    parse_dictionary["maxthreads"] = lambda parser: parser.add_argument(
        "-n_thread_fam",
        default=1,
        required=False,
        type=int,
        help="Maximum number of family threads to use. Default: 1.",
    )
    parse_dictionary["n_thread_loci"] = lambda parser: parser.add_argument(
        "-n_thread_loci",
        default=1,
        required=False,
        type=_positive_int,
        help=("Number of locus threads to use inside each family peel. Default: 1."),
    )
    return parse_dictionary


def _add_program_arguments(parser):
    """Add general program options."""

    program_parser = parser.add_argument_group("Program options")
    program_parser.add_argument(
        "-h",
        "-help",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    program_parser.add_argument(
        "-version",
        default=None,
        action="version",
        version="%(prog)s " + version_version,
        help="Show program's version number and exit.",
    )


def _add_input_individual_arguments(parser):
    """Add input options for individual-level data."""

    input_parser = parser.add_argument_group("Input options: individuals")
    add_arguments_from_dictionary(
        input_parser,
        get_input_options(),
        options=[
            "pedigree",
            "genotypes",
            "seqfile",
            "phasefile",
            "phenotype",
        ],
    )
    input_parser.add_argument(
        "-x_chr",
        action="store_true",
        required=False,
        help="Indicate that input data is for the X chromosome (see details in the docs).",
    )
    input_parser.add_argument(
        "-phased_geno_prob_file",
        default=None,
        required=False,
        type=str,
        nargs="*",
        help="Optional external phased genotype probability file(s) "
        "(see format details in the docs). "
        "This will provide the starting internal genotype probability state. ",
    )


def _add_input_marker_arguments(parser):
    """Add input options for marker ranges and maps."""

    marker_parser = parser.add_argument_group("Input options: markers")
    marker_parser.add_argument(
        "-map_file",
        default=None,
        required=False,
        type=str,
        help="Map file for loci in genomic data files (see format details in the docs).",
    )
    add_arguments_from_dictionary(
        marker_parser,
        get_input_options(),
        options=["startsnp", "stopsnp"],
    )


def _add_input_model_parameter_arguments(parser):
    """Add input options for model parameters."""

    parameter_parser = parser.add_argument_group(
        "Input options: model parameters and other"
    )
    add_arguments_from_dictionary(
        parameter_parser,
        get_input_options(),
        options=[
            "alt_allele_prob_file",
            "main_metafounder",
            "pheno_penetrance_prob_file",
        ],
    )
    parameter_parser.add_argument(
        "-rec_length",
        default=1.0,
        required=False,
        type=float,
        help="Recombination length of the chromosome in Morgans. Default: 1.00",
    )
    add_arguments_from_dictionary(
        parameter_parser,
        get_probability_options(),
        options=["mut_prob", "geno_error_prob", "seq_error_prob"],
    )


def _add_output_individual_arguments(parser):
    """Add output options for individual-level result files."""

    output_parser = parser.add_argument_group("Output options: individuals")
    output_parser.add_argument(
        "-no_dosage",
        action="store_true",
        required=False,
        help="Suppress default output of allele dosages (see format details in the docs).",
    )
    output_parser.add_argument(
        "-geno",
        action="store_true",
        required=False,
        help="Call and output genotypes (see format details in the docs). "
        "The default genotype calling threshold is set to 1/3.",
    )
    output_parser.add_argument(
        "-geno_threshold",
        default=None,
        required=False,
        type=float,
        nargs="*",
        help="Custom genotype calling threshold(s) from the genotype probabilities. "
        "Multiple space separated values allowed.\
        Value(s) less than 1/3 are replaced by 1/3.",
    )
    output_parser.add_argument(
        "-geno_prob",
        action="store_true",
        required=False,
        help="Output genotype probabilities (see format details in the docs).",
    )
    output_parser.add_argument(
        "-phased_geno_prob",
        action="store_true",
        required=False,
        help="Output phased genotype probabilities (see format details in the docs).",
    )
    output_parser.add_argument(
        "-hap",
        action="store_true",
        required=False,
        help="Call and output haplotypes (see format details in the docs). "
        "The default haplotype calling threshold is set to 1/2.",
    )
    output_parser.add_argument(
        "-hap_threshold",
        default=None,
        required=False,
        type=float,
        nargs="*",
        help="Custom haplotype calling threshold(s) from the phased genotype probabilities. "
        "Multiple space separated values allowed.\
        Value(s) less than 1/2 are replaced by 1/2.",
    )
    output_parser.add_argument(
        "-seg_prob",
        action="store_true",
        required=False,
        help="Output segregation probabilities (see format details in the docs).",
    )
    output_parser.add_argument(
        "-pheno_prob",
        action="store_true",
        required=False,
        help="Output phenotype probabilities (see format details in the docs).",
    )
    output_parser.add_argument(
        "-rec_prob",
        action="store_true",
        required=False,
        help="Output recombination probabilities (NOT IMPLEMENTED YET).",
    )
    output_parser.add_argument(
        "-alt_allele_prob",
        action="store_true",
        required=False,
        help="Output alternative allele frequencies (see format details in the docs). "
        "Output 0.5 if none of ``est_start_alt_allele_prob``, ``est_alt_allele_prob``, "
        "or ``alt_allele_prob_file`` is used.",
    )
    output_parser.add_argument(
        "-pheno_penetrance_prob",
        action="store_true",
        required=False,
        help="Output phenotype penetrance probabilities (see format details in the docs).",
    )


def _add_output_io_arguments(parser):
    """Add output prefix, ordering, and formatting options."""

    output_parser = parser.add_argument_group("Output options: prefix, order, and IO")
    output_parser.add_argument(
        "-out_file",
        required=True,
        type=str,
        help='The output file prefix. All file outputs will be named as "PREFIX.OUTPUT.txt", '
        'where "OUTPUT" is the type of output (for example, "dosage" and "geno_prob").',
    )
    add_arguments_from_dictionary(
        output_parser,
        get_output_options(),
        options=["writekey", "onlykeyed", "out_digits"],
    )


def _add_peeling_method_arguments(parser):
    """Add peeling strategy and hybrid second-stage options."""

    method_parser = parser.add_argument_group("Peeling methods: strategy")
    method_parser.add_argument(
        "-method",
        default=None,
        required=False,
        type=str,
        help="Peeling method: single or multi. Default: multi.",
    )

    hybrid_parser = parser.add_argument_group("Peeling methods: hybrid second stage")
    hybrid_parser.add_argument(
        "-seg_file",
        default=None,
        required=False,
        type=str,
        help="Segregation probabilities file (see format details in the docs).",
    )
    hybrid_parser.add_argument(
        "-seg_map_file",
        default=None,
        required=False,
        type=str,
        help="Map file for loci in the segregation probabilities file.",
    )


def _add_peeling_parameter_arguments(parser):
    """Add peeling runtime and estimation options."""

    computational_parser = parser.add_argument_group(
        "Peeling parameters: computational"
    )
    computational_parser.add_argument(
        "-n_cycle",
        default=5,
        required=False,
        type=int,
        help="Number of peeling cycles. Default: 5.",
    )
    add_arguments_from_dictionary(
        computational_parser,
        get_multithread_options(),
        options=["maxthreads", "n_thread_loci"],
    )

    estimation_parser = parser.add_argument_group(
        "Peeling parameters: model estimation"
    )
    estimation_parser.add_argument(
        "-est_start_alt_allele_prob",
        action="store_true",
        required=False,
        help="Estimate from all inputted genomic data prior to peeling "
        "and output alternative allele probabilities "
        "(see format details in the docs).",
    )
    estimation_parser.add_argument(
        "-est_alt_allele_prob",
        action="store_true",
        required=False,
        help="Estimate after each peeling cycle and output alternative allele probabilities "
        "(see format details in the docs).",
    )
    estimation_parser.add_argument(
        "-est_geno_error_prob",
        action="store_true",
        required=False,
        help="Estimate after each peeling cycle and output genotype error probabilities "
        "(see format details in the docs).",
    )
    estimation_parser.add_argument(
        "-est_seq_error_prob",
        action="store_true",
        required=False,
        help="Estimate after each peeling cycle and output sequence error probabilities "
        "(see format details in the docs).",
    )
    estimation_parser.add_argument(
        "-est_pheno_penetrance_prob",
        action="store_true",
        required=False,
        help="Estimate after each peeling cycle and output phenotype penetrance probabilities "
        "(see format details in the docs).",
    )
    estimation_parser.add_argument(
        "-no_phase_founder",
        action="store_true",
        required=False,
        help="Suppress phasing a heterozygous allele (if such an allele can be found) "
        "in genotyped individuals without genotyped parents.",
    )


def build_arg_parser():
    """Build and return the AlphaPeel argument parser.

    :return: AlphaPeel argument parser.
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(description="", add_help=False)

    _add_program_arguments(parser)
    _add_input_individual_arguments(parser)
    _add_input_marker_arguments(parser)
    _add_input_model_parameter_arguments(parser)
    _add_output_individual_arguments(parser)
    _add_output_io_arguments(parser)
    _add_peeling_method_arguments(parser)
    _add_peeling_parameter_arguments(parser)

    return parser


def _argv_contains_version(argv):
    """Return whether argv asks only for the AlphaPeel version."""

    args = sys.argv[1:] if argv is None else list(argv)
    return "-version" in args


def _zero_based_locus(value):
    """Convert a user-facing 1-based locus value to Python indexing."""

    if value is None:
        return None

    return value - 1


def _positive_int(value):
    """Parse a positive integer command-line argument."""

    parsed_value = int(value)
    if parsed_value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed_value


def normalise_alphapeel_args(args):
    """Populate tinyhouse-compatible argument names for AlphaPeel."""

    for tinyhouse_name, alphapeel_name in ALPHAPEEL_ARGUMENT_ALIASES.items():
        setattr(args, tinyhouse_name, getattr(args, alphapeel_name))

    args.startsnp = _zero_based_locus(args.start_snp)
    args.stopsnp = _zero_based_locus(args.stop_snp)

    return args


def parse_alphapeel_args(parser, argv=None):
    """Parse AlphaPeel arguments and normalise names used by tinyhouse."""

    if _argv_contains_version(argv):
        parser.parse_args(sys.argv[1:] if argv is None else list(argv))
        sys.exit(0)

    args = parseArgs("AlphaPeel", parser, argv=argv)

    return normalise_alphapeel_args(args)


def get_args(argv=None):
    """Presents and collects the arguments from the command line.

    :return: the user input arguments for the AlphaPeel program
    :rtype: argparse.Namespace
    """

    return parse_alphapeel_args(build_arg_parser(), argv=argv)


def main(argv=None):
    """Main function for the AlphaPeel program.
    This function collects the arguments from the command line and runs the peeling algorithm.
    """
    docs_link = f"https://alphapeel.readthedocs.io/en/v{version_version}/usage.html"
    print_boilerplate("AlphaPeel", version=version_version, docs=docs_link)
    args = get_args(argv=argv)

    pedigree = Pedigree()
    readInPedigreeFromInputs(pedigree, args)

    single_locus_mode = args.method == "single"
    if args.method == "multi" and args.segfile:
        warnings.warn("Running in multi-locus mode, external segfile ignored")

    # For now, only support a single phenotype (will extend in future)
    if args.phenotype is not None and pedigree.nPheno > 1:
        warnings.warn(
            "Currently only a single phenotype is supported. Phenotype information will be ignored."
        )
        pedigree.phenoPenetrance = None
        args.phenoPenetrance = None
        for ind in pedigree:
            ind.phenotype = None

    peeling_info = create_peeling_info(
        pedigree, args, phase_founder=(not args.no_phase_founder)
    )

    if single_locus_mode:
        print("Generating seg estimates")
        generate_single_locus_segregation(peeling_info, pedigree, args)
    run_peeling_cycles(
        pedigree, peeling_info, args, single_locus_mode=single_locus_mode
    )

    write_requested_outputs(
        pedigree, peeling_info, args, single_locus_mode=single_locus_mode
    )


if __name__ == "__main__":
    main()
