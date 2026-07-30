"""Peeling information setup and the JIT container used by peeling cycles."""

import warnings
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from numba import jit, optional, boolean, int8, uint32, float32
from numba.experimental import jitclass

from .peeling_io import add_penetrance_from_external_file
from ..tinyhouse import InputOutput
from ..tinyhouse.ProbMath import (
    generateSegregation,
    generateSegregationXYChrom,
    generateSegregationXXChrom,
    getGenotypeProbabilities,
    updateGenoProbsFromPhenotype,
)
from ..tinyhouse.HaplotypeOperations import ind_fillInGenotypesFromPhase


@dataclass
class PeelingCycleContext:
    """Prepared state reused by one or more peeling cycles."""

    pedigree: object
    peeling_info: object
    n_fam_threads: int
    single_locus_mode: bool
    jit_generations: list
    locus_thread_blocks: list


def create_peeling_info(pedigree, args, phase_founder=False):
    """Create the peeling information object and initialize model probabilities.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :param phase_founder: whether to phase genotyped founders using heterozygous loci,
        defaults to False
    :type phase_founder: bool, optional
    :return: peeling_info: a peeling information object containing all the
        necessary information for peeling
    :rtype: JitPeelingInformation
    """
    n_loci = pedigree.nLoci

    peeling_info = JitPeelingInformation(
        n_ind=pedigree.maxIdn, n_fam=pedigree.maxFam, n_loci=n_loci
    )

    initialize_peeling_info_model(peeling_info, args)
    initialize_individual_probabilities(pedigree, peeling_info, args, phase_founder)
    add_external_phased_genotype_probabilities(pedigree, peeling_info, args)

    return peeling_info


def create_locus_block_peeling_infos(peeling_info, n_blocks):
    """Create view-backed peeling-info objects for contiguous locus blocks."""

    blocks = []
    n_loci = peeling_info.n_loci
    n_blocks = min(n_blocks, n_loci)
    chunk_size = (n_loci + n_blocks - 1) // n_blocks

    for start in range(0, n_loci, chunk_size):
        stop = min(start + chunk_size, n_loci)
        blocks.append(
            (start, stop, view_peeling_info_locus_block(peeling_info, start, stop))
        )

    return blocks


def view_peeling_info_locus_block(peeling_info, start, stop):
    """Create a peeling-info object whose arrays are views of one locus block."""

    block_n_loci = stop - start
    # Seed the jitclass with the smallest valid allocation, then replace fields
    # with views. This avoids allocating full block-sized arrays just to discard
    # them immediately.
    block_info = JitPeelingInformation(n_ind=1, n_fam=1, n_loci=1)

    block_info.n_ind = peeling_info.n_ind
    block_info.n_fam = peeling_info.n_fam
    block_info.n_loci = block_n_loci
    block_info.iteration = peeling_info.iteration
    block_info.is_x_chr = peeling_info.is_x_chr
    # The fields below are declared on JitPeelingInformation. For block views,
    # we intentionally rebind them from allocated arrays to slices of the full
    # peeling_info arrays.
    # pylint: disable=attribute-defined-outside-init
    block_info.sex = peeling_info.sex

    block_info.anterior = peeling_info.anterior[:, :, start:stop]
    block_info.posterior = peeling_info.posterior[:, :, start:stop]
    block_info.penetrance = peeling_info.penetrance[:, :, start:stop]
    block_info.segregation = peeling_info.segregation[:, :, start:stop]
    block_info.posterior_sire_contribution = peeling_info.posterior_sire_contribution[
        :, :, start:stop
    ]
    block_info.posterior_dam_contribution = peeling_info.posterior_dam_contribution[
        :, :, start:stop
    ]

    block_info.geno_error = peeling_info.geno_error[start:stop]
    block_info.seq_error = peeling_info.seq_error[start:stop]
    view_block_transmission(peeling_info, block_info, start, stop)
    view_block_positions(peeling_info, block_info, start, stop)
    copy_shared_segregation_tensors(peeling_info, block_info)

    return block_info


def view_block_transmission(peeling_info, block_info, start, stop):
    """Point a block at the full transmission rates inside the block."""

    if block_info.n_loci <= 1:
        return

    block_info.transmission_rate = peeling_info.transmission_rate[start : stop - 1]


def view_block_positions(peeling_info, block_info, start, stop):
    """Point a block at locus positions when present."""

    block_info.positions = None
    if peeling_info.positions is not None:
        block_info.positions = peeling_info.positions[start:stop]


def copy_shared_segregation_tensors(peeling_info, block_info):
    """Copy chromosome-level segregation tensor references to a block."""

    block_info.segregation_tensor = peeling_info.segregation_tensor
    block_info.segregation_tensor_norm = peeling_info.segregation_tensor_norm
    block_info.segregation_tensor_xy = peeling_info.segregation_tensor_xy
    block_info.segregation_tensor_xy_norm = peeling_info.segregation_tensor_xy_norm
    block_info.segregation_tensor_xx = peeling_info.segregation_tensor_xx
    block_info.segregation_tensor_xx_norm = peeling_info.segregation_tensor_xx_norm


def initialize_peeling_info_model(peeling_info, args):
    """Initialize chromosome, map, error, transmission, and segregation settings."""

    peeling_info.is_x_chr = args.x_chr
    set_map_positions(peeling_info, args)
    setup_segregation_tensors(peeling_info, args.mut_prob)
    peeling_info.geno_error[:] = args.geno_error_prob
    peeling_info.seq_error[:] = args.seq_error_prob
    setup_transmission(args.rec_length, peeling_info)


def set_map_positions(peeling_info, args):
    """Set marker positions when map-aware transmission rates are used."""

    # Positions are only needed when map-aware transmission rates are used.
    peeling_info.positions = None
    if not args.map_file:
        return

    if args.stopsnp is not None:
        positions = InputOutput.readMapFile(
            args.map_file, args.startsnp, args.stopsnp + 1
        )[2]
    else:
        positions = InputOutput.readMapFile(args.map_file)[2]
    peeling_info.positions = np.array(positions, dtype=np.uint32)


def setup_segregation_tensors(peeling_info, mut_prob):
    """Set segregation tensors used by the peeling kernels."""

    # Segregation tensors encode P(parent genotypes, child genotype, segregation state).
    # The *_norm tensors are partial tensors used as normalizing constants.
    peeling_info.segregation_tensor = generateSegregation(mu=mut_prob)
    peeling_info.segregation_tensor_norm = generateSegregation(
        mu=mut_prob, partial=True
    )
    if peeling_info.is_x_chr:
        peeling_info.segregation_tensor_xy = generateSegregationXYChrom(mu=mut_prob)
        peeling_info.segregation_tensor_xy_norm = generateSegregationXYChrom(
            mu=mut_prob, partial=True
        )
        peeling_info.segregation_tensor_xx = generateSegregationXXChrom(mu=mut_prob)
        peeling_info.segregation_tensor_xx_norm = generateSegregationXXChrom(
            mu=mut_prob, partial=True
        )


def initialize_individual_probabilities(pedigree, peeling_info, args, phase_founder):
    """Initialize per-individual sex, penetrance, and X-chromosome segregation."""

    for ind in pedigree:
        peeling_info.sex[ind.idn] = ind.sex

        fill_genotypes_from_phase(ind)
        ind_penetrance = get_individual_penetrance(
            ind, pedigree, peeling_info, args, phase_founder
        )
        initialize_x_chr_segregation(ind, peeling_info)
        peeling_info.penetrance[ind.idn, :, :] = ind_penetrance


def fill_genotypes_from_phase(ind):
    """Fill genotypes from phase when both genotype and haplotype inputs exist."""

    if ind.genotypes is not None and ind.haplotypes is not None:
        ind_fillInGenotypesFromPhase(ind)


def get_individual_penetrance(ind, pedigree, peeling_info, args, phase_founder):
    """Build one individual's initial genotype penetrance matrix."""

    ind_penetrance = get_base_individual_penetrance(ind, peeling_info)
    ind_penetrance = apply_phenotype_to_penetrance(ind, pedigree, ind_penetrance)
    apply_founder_phase_penetrance(
        ind, peeling_info, args, phase_founder, ind_penetrance
    )
    return ind_penetrance


def get_base_individual_penetrance(ind, peeling_info):
    """Build genotype/read penetrance before phenotype or founder-phase adjustments."""

    x_chr_male_flag = peeling_info.is_x_chr and ind.sex == 0
    return getGenotypeProbabilities(
        peeling_info.n_loci,
        ind.genotypes,
        ind.reads,
        peeling_info.geno_error,
        peeling_info.seq_error,
        x_chr_male_flag,
    )


def initialize_x_chr_segregation(ind, peeling_info):
    """Set initial X-chromosome segregation probabilities for one individual."""

    if not peeling_info.is_x_chr:
        return

    ind_segregation = peeling_info.segregation[ind.idn, :, :]
    if ind.sex == 0:
        # Males use the paternal X states.
        ind_segregation[0, :] = 0.5
        ind_segregation[1, :] = 0.5
        ind_segregation[2, :] = 0
        ind_segregation[3, :] = 0
    elif ind.sex == 1:
        # Females use the maternal X states.
        ind_segregation[0, :] = 0
        ind_segregation[1, :] = 0
        ind_segregation[2, :] = 0.5
        ind_segregation[3, :] = 0.5


def apply_phenotype_to_penetrance(ind, pedigree, ind_penetrance):
    """Apply phenotype probabilities to one individual's penetrance."""

    if ind.phenotype is None:
        return ind_penetrance

    # Phenotypes further weight genotype penetrance.
    # Phenotype updates currently assume the existing single-locus phenotype model;
    # multi-phenotype or multi-locus phenotype handling needs an explicit extension.
    return updateGenoProbsFromPhenotype(
        ind_penetrance,
        ind.phenotype,
        pedigree.phenoPenetrance,
    )


def apply_founder_phase_penetrance(
    ind, peeling_info, args, phase_founder, ind_penetrance
):
    """Set phased-founder penetrance at the heterozygous midpoint when requested."""

    if not (ind.isGenotypedFounder() and phase_founder and ind.genotypes is not None):
        return

    loci = get_het_midpoint(ind.genotypes)
    if loci is None or (peeling_info.is_x_chr and ind.sex == 0):
        return

    error = args.geno_error_prob
    ind_penetrance[:, loci] = np.array(
        [error / 3, error / 3, 1 - error, error / 3], dtype=np.float32
    )


def add_external_phased_genotype_probabilities(pedigree, peeling_info, args):
    """Apply external phased genotype probability files when configured."""

    if args.phased_geno_prob_file is None:
        return

    warn_for_external_phased_genotype_options(peeling_info, args)
    for pen in args.phased_geno_prob_file:
        add_penetrance_from_external_file(pedigree, peeling_info, pen, args)


def warn_for_external_phased_genotype_options(peeling_info, args):
    """Warn and disable incompatible options for external phased probabilities."""

    if peeling_info.is_x_chr:
        warnings.warn(
            "Using an external phased genotype probability file and "
            "the x_chr option is highly discouraged. Please do not use.",
            UserWarning,
        )

    if args.est_geno_error_prob:
        warnings.warn(
            "External phased genotype probability file included, "
            "but est_geno_error_prob flag used. "
            "The two options are incompatible. est_geno_error_prob set to false.",
            UserWarning,
        )
        args.est_geno_error_prob = False

    if args.est_seq_error_prob:
        warnings.warn(
            "External phased genotype probability file included, "
            "but est_seq_error_prob flag used. "
            "The two options are incompatible. est_seq_error_prob set to false.",
            UserWarning,
        )
        args.est_seq_error_prob = False


def setup_transmission(length, peeling_info):
    """Set transmission rates from the distance between neighbouring loci.

    If no map positions are provided, loci are spaced evenly along the chromosome.
    Otherwise positions are rescaled to [0, 1] and multiplied by chromosome length.

    :param length: Estimated recombination length of the chromosome in Morgans. Default = 1.00
    :type length: float
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`peeling_info_module.JitPeelingInformation`
    :return: None. Updates peeling_info.transmission_rate in place
    """
    # Relative positions of the loci on a chromosome, scaled between 0 and 1.
    if peeling_info.positions is None:
        local_map = np.linspace(0, 1, num=peeling_info.n_loci, dtype=np.float32)
    else:
        local_map = (peeling_info.positions - peeling_info.positions[0]) / (
            peeling_info.positions[-1] - peeling_info.positions[0]
        )
    for i in range(peeling_info.n_loci - 1):
        distance = local_map[i + 1] - local_map[i]
        distance = distance * length
        peeling_info.transmission_rate[i] = distance


@jit(nopython=True)
def get_het_midpoint(geno):
    """Finds the midpoint of the heterozygous loci in a genotype array.

    :param geno: observed genotypes for an individual collected via user input.
    :type geno: 1D numpy array of Int8 with length n_loci
    :return: The index of the first heterozygous locus found,
        or None if no heterozygous loci are present.
    :rtype: int or None
    """
    n_loci = len(geno)
    midpoint = int(n_loci / 2)
    index = 0
    while index < n_loci / 2:
        if midpoint + index < n_loci:
            if geno[midpoint + index] == 1:
                return midpoint + index
        if midpoint - index >= 0:
            if geno[midpoint - index] == 1:
                return midpoint - index
        index += 1
    return None


spec = OrderedDict()
spec["n_ind"] = uint32
spec["n_fam"] = uint32
spec["n_loci"] = uint32

spec["is_x_chr"] = boolean
spec["sex"] = int8[:]

# Individual terms: Each will be n_ind x 4 x n_loci
spec["anterior"] = float32[:, :, :]
spec["posterior"] = float32[:, :, :]
spec["penetrance"] = float32[:, :, :]
spec["segregation"] = optional(float32[:, :, :])

# Family posterior contributions. Each will be n_fam x 4 x n_loci.
spec["posterior_sire_contribution"] = float32[:, :, :]
spec["posterior_dam_contribution"] = float32[:, :, :]

# Segregation tensors. Each of these will be either 4x4x4x4 or 4x4x4
spec["segregation_tensor"] = optional(float32[:, :, :, :])
spec["segregation_tensor_norm"] = optional(float32[:, :, :])
spec["segregation_tensor_xx"] = optional(float32[:, :, :, :])
spec["segregation_tensor_xy"] = optional(float32[:, :, :, :])
spec["segregation_tensor_xx_norm"] = optional(float32[:, :, :])
spec["segregation_tensor_xy_norm"] = optional(float32[:, :, :])

# Marker specific rates:
spec["geno_error"] = optional(float32[:])
spec["seq_error"] = optional(float32[:])
spec["transmission_rate"] = optional(float32[:])

spec["positions"] = optional(uint32[:])
spec["iteration"] = uint32


@jitclass(spec)
# This jitclass is the compiled peeling state container; keeping the arrays as
# direct attributes preserves clear numba field types and call-site access.
# pylint: disable=too-many-instance-attributes
class JitPeelingInformation:
    """Holds the peeling information for a given pedigree.


    :param object: peeling information object
    :type object: class:`JitPeelingInformation`
    """

    def __init__(self, n_ind, n_fam, n_loci):
        """Initialize the peeling information object.

        :param n_ind: number of individuals in the pedigree
        :type n_ind: int
        :param n_fam: number of families in the pedigree
        :type n_fam: int
        :param n_loci: number of loci in genotype input
        :type n_loci: int
        """
        self.iteration = 0
        self.n_ind = n_ind
        self.n_fam = n_fam
        self.n_loci = n_loci

        self.is_x_chr = False

        self.construct()

        # These are filled in from create_peeling_info, above.
        self.positions = None
        self.segregation_tensor = None
        self.segregation_tensor_norm = None

        self.segregation_tensor_xy = None
        self.segregation_tensor_xy_norm = None
        self.segregation_tensor_xx = None
        self.segregation_tensor_xx_norm = None

    def construct(self):
        """Allocate the probability arrays with uniform starting values."""
        base_value = 0.25
        self.sex = np.full(self.n_ind, 0, dtype=np.int8)

        self.anterior = np.full(
            (self.n_ind, 4, self.n_loci), base_value, dtype=np.float32
        )
        self.posterior = np.full(
            (self.n_ind, 4, self.n_loci), base_value, dtype=np.float32
        )
        self.penetrance = np.full(
            (self.n_ind, 4, self.n_loci), base_value, dtype=np.float32
        )

        self.segregation = np.full(
            (self.n_ind, 4, self.n_loci), base_value, dtype=np.float32
        )

        self.posterior_sire_contribution = np.full(
            (self.n_fam, 4, self.n_loci), base_value, dtype=np.float32
        )
        self.posterior_dam_contribution = np.full(
            (self.n_fam, 4, self.n_loci), base_value, dtype=np.float32
        )

        self.geno_error = np.full((self.n_loci), 0, dtype=np.float32)
        self.seq_error = np.full((self.n_loci), 0, dtype=np.float32)
        self.transmission_rate = np.full((self.n_loci - 1), 0, dtype=np.float32)

    def get_geno_probs(self, idn, sex=None):
        """Estimates the genotype probabilities for a given individual.

        :param idn: Internal number for an individual in the pedigree.
        :type idn: int
        :param sex: 0 is male; 1 is female
        :type sex: int
        :return: geno_probs: the genotype probabilities for the individual
        :rtype: 2D numpy array of float32 with shape 4 x n_loci
        """
        anterior = self.anterior[idn, :, :]
        posterior = self.posterior[idn, :, :]
        penetrance = self.penetrance[idn, :, :]
        geno_probs = anterior * posterior * penetrance
        if self.is_x_chr and sex == 0:  # male
            geno_probs0 = geno_probs[0, :]
            geno_probs1 = geno_probs[1, :]
            geno_probs2 = geno_probs[2, :]
            geno_probs3 = geno_probs[3, :]
            geno_probs0 += geno_probs2
            geno_probs3 += geno_probs1
            geno_probs1[:] = 0
            geno_probs2[:] = 0

        geno_probs /= np.sum(geno_probs, 0)
        return geno_probs

    def get_pheno_probs(self, idn, pheno_penetrance):
        """Estimates the phenotype probabilities for a given individual.

        :param idn: Internal number for an individual in the pedigree.
        :type idn: int
        :param pheno_penetrance: phenotype penetrance values
        :type pheno_penetrance: 2D numpy array of float32 with shape
            4 x number of phenotype categories
        :return: pheno_probs: the phenotype probabilities for the individual
        :rtype: 2D numpy array of float32 with shape number of phenotype categories x 1
        """
        geno_probs = self.get_geno_probs(idn)
        rg_pheno = pheno_penetrance.shape[1]
        i = 0
        pheno_probs = np.zeros((rg_pheno, 1), dtype=np.float32)
        while i < rg_pheno:
            total = 0.0
            for genotype in range(4):
                penetrance = pheno_penetrance[genotype, i]
                for locus in range(self.n_loci):
                    total += geno_probs[genotype, locus] * penetrance
            pheno_probs[i, 0] = total
            i += 1

        pheno_probs /= np.sum(pheno_probs, 0)
        return pheno_probs
