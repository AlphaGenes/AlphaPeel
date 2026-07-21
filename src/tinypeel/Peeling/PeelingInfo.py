"""Peeling information setup and the JIT container used by peeling cycles."""

from numba import jit, optional, boolean, int8, uint32, float32
from numba.experimental import jitclass
import numpy as np
from collections import OrderedDict
import warnings

from ..tinyhouse import ProbMath
from ..tinyhouse import HaplotypeOperations
from ..tinyhouse import InputOutput


def create_peeling_info(pedigree, args, phase_founder=False):
    """Create the peeling information object and initialize model probabilities.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param args: argument container with configuration options for peeling
    :type args: argparse.Namespace or similar object with attributes
    :param phase_founder: whether to phase genotyped founders using heterozygous loci, defaults to False
    :type phase_founder: bool, optional
    :return: peeling_info: a peeling information object containing all the necessary information for peeling
    :rtype: jit_peeling_information
    """
    n_loci = pedigree.nLoci

    peeling_info = jit_peeling_information(
        n_ind=pedigree.maxIdn, n_fam=pedigree.maxFam, n_loci=n_loci
    )

    peeling_info.is_x_chr = args.x_chr
    # Positions are only needed when map-aware transmission rates are used.
    peeling_info.positions = None
    if args.map_file:
        if args.stopsnp is not None:
            peeling_info.positions = np.array(
                InputOutput.readMapFile(args.map_file, args.startsnp, args.stopsnp + 1)[
                    2
                ],
                dtype=np.uint32,
            )
        else:
            peeling_info.positions = np.array(
                InputOutput.readMapFile(args.map_file)[2], dtype=np.uint32
            )

    mut_prob = args.mut_prob

    # Segregation tensors encode P(parent genotypes, child genotype, segregation state).
    # The *_norm tensors are partial tensors used as normalizing constants.
    peeling_info.segregation_tensor = ProbMath.generateSegregation(mu=mut_prob)
    peeling_info.segregation_tensor_norm = ProbMath.generateSegregation(
        mu=mut_prob, partial=True
    )
    if peeling_info.is_x_chr:
        peeling_info.segregation_tensor_xy = ProbMath.generateSegregationXYChrom(
            mu=mut_prob
        )
        peeling_info.segregation_tensor_xy_norm = ProbMath.generateSegregationXYChrom(
            mu=mut_prob, partial=True
        )
        peeling_info.segregation_tensor_xx = ProbMath.generateSegregationXXChrom(
            mu=mut_prob
        )
        peeling_info.segregation_tensor_xx_norm = ProbMath.generateSegregationXXChrom(
            mu=mut_prob, partial=True
        )

    peeling_info.geno_error[:] = args.geno_error_prob
    peeling_info.seq_error[:] = args.seq_error_prob
    setup_transmission(args.rec_length, peeling_info)

    for ind in pedigree:
        peeling_info.sex[ind.idn] = ind.sex

        if ind.genotypes is not None and ind.haplotypes is not None:
            HaplotypeOperations.ind_fillInGenotypesFromPhase(ind)

        x_chr_male_flag = peeling_info.is_x_chr and ind.sex == 0

        ind_penetrance = ProbMath.getGenotypeProbabilities(
            peeling_info.n_loci,
            ind.genotypes,
            ind.reads,
            peeling_info.geno_error,
            peeling_info.seq_error,
            x_chr_male_flag,
        )

        if peeling_info.is_x_chr:
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
        if ind.phenotype is not None:
            # Phenotypes further weight genotype penetrance.
            # Phenotype updates currently assume the existing single-locus phenotype model;
            # multi-phenotype or multi-locus phenotype handling needs an explicit extension.
            ind_penetrance = ProbMath.updateGenoProbsFromPhenotype(
                ind_penetrance,
                ind.phenotype,
                pedigree.phenoPenetrance,
            )

        if ind.isGenotypedFounder() and phase_founder and ind.genotypes is not None:
            loci = get_het_midpoint(ind.genotypes)
            if loci is not None:
                error = args.geno_error_prob
                if (not peeling_info.is_x_chr) or (ind.sex != 0):
                    ind_penetrance[:, loci] = np.array(
                        [error / 3, error / 3, 1 - error, error / 3], dtype=np.float32
                    )

        peeling_info.penetrance[ind.idn, :, :] = ind_penetrance

    if args.phased_geno_prob_file is not None:
        if peeling_info.is_x_chr:
            warnings.warn(
                "Using an external phased genotype probability file and the x_chr option is highly discouraged. Please do not use.",
                UserWarning,
            )

        if args.est_geno_error_prob:
            warnings.warn(
                "External phased genotype probability file included, but est_geno_error_prob flag used. The two options are incompatible. est_geno_error_prob set to false.",
                UserWarning,
            )
            args.est_geno_error_prob = False

        if args.est_seq_error_prob:
            warnings.warn(
                "External phased genotype probability file included, but est_seq_error_prob flag used. The two options are incompatible. est_seq_error_prob set to false.",
                UserWarning,
            )
            args.est_seq_error_prob = False

        for pen in args.phased_geno_prob_file:
            add_penetrance_from_external_file(pedigree, peeling_info, pen, args)

    return peeling_info


def setup_transmission(length, peeling_info):
    """Set transmission rates from the distance between neighbouring loci.

    If no map positions are provided, loci are spaced evenly along the chromosome.
    Otherwise positions are rescaled to [0, 1] and multiplied by chromosome length.

    :param length: Estimated recombination length of the chromosome in Morgans. Default = 1.00
    :type length: float
    :param peeling_info: Peeling information container.
    :type peeling_info: class:`PeelingInfo.jit_peeling_information`
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


def add_penetrance_from_external_file(pedigree, peeling_info, file_name, args):
    """Read external genotype penetrance values and multiply them into penetrance.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param peeling_info: Peeling information container
    :type peeling_info: class:`PeelingInfo.jit_peeling_information`
    :param file_name: path to the external penetrance file
    :type file_name: str
    :param args: argument container with configuration options for peeling set up, including startsnp and stopsnp.
    :type args: argparse.Namespace or similar object with attributes
    :return: None. Updates peeling_info.penetrance in place
    """
    print("Reading in penetrance file:", file_name)
    with open(file_name) as f:
        e = 0
        for line in f:
            parts = line.split()
            idx = parts[0]
            parts = parts[1:]

            if args.startsnp is not None:
                parts = parts[args.startsnp : args.stopsnp + 1]

            penetrance_line = np.array([float(val) for val in parts], dtype=np.float32)

            if idx not in pedigree.individuals:
                warnings.warn(
                    "Individual",
                    idx,
                    "not found in pedigree. Individual ignored.",
                    UserWarning,
                )
            else:
                ind = pedigree.individuals[idx]
                penetrance_prob = peeling_info.penetrance[ind.idn, e, :]
                penetrance_prob *= penetrance_line
                e = (e + 1) % 4


@jit(nopython=True)
def get_het_midpoint(geno):
    """Finds the midpoint of the heterozygous loci in a genotype array.

    :param geno: observed genotypes for an individual collected via user input.
    :type geno: 1D numpy array of Int8 with length n_loci
    :return: The index of the first heterozygous locus found, or None if no heterozygous loci are present.
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
class jit_peeling_information(object):
    """Holds the peeling information for a given pedigree.


    :param object: peeling information object
    :type object: class:`jit_peeling_information`
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
        :type pheno_penetrance: 2D numpy array of float32 with shape 4 x number of phenotype categories
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
