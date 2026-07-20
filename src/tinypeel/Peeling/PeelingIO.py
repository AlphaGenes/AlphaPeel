import numpy as np
from numba import jit
from contextlib import ExitStack
from ..tinyhouse import InputOutput


def write_out_parameters(peeling_info):
    """Writes out the geno error rate, seq error rate, and recombination probabilities.

    :param peeling_info: Peeling information container
    :type peeling_info: class:`PeelingInfo.jit_peeling_information`
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    if args.est_geno_error_prob:
        np.savetxt(
            args.out_file + ".geno_error_prob.txt",
            peeling_info.geno_error,
        )
    if args.est_seq_error_prob:
        np.savetxt(
            args.out_file + ".seq_error_prob.txt",
            peeling_info.seq_error,
        )
    if args.rec_prob:
        np.savetxt(
            args.out_file + ".rec_prob.txt",
            np.empty((1, 1)),
        )  # not implemented atm, just as a placeholder


def write_out_alt_allele_prob(pedigree):
    """Writes out the alternative allele probabilities for each locus and metafounder in the pedigree.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :return: None. Writes to a file specified in the InputOutput.args.
    """
    args = InputOutput.args

    # Custom sorting key to extract numeric part of metafounder keys if it's an integer
    def sort_key(mf_key):
        part = mf_key.split("_")[1]
        return (0, int(part)) if part.isdigit() else (1, part)

    # Order the metafounders: MF_1, MF_2, ..., MF_11, etc., otherwise keep original order
    sorted_aap = dict(sorted(pedigree.AAP.items(), key=lambda item: sort_key(item[0])))
    sorted_mf = list(sorted_aap.keys())
    # Combine data into a single 2D array
    combined_aap = np.hstack(
        [sorted_aap[key].reshape(pedigree.nLoci, -1) for key in sorted_mf]
    )
    # Save into text file with metafounders heading columns
    np.savetxt(
        args.out_file + ".alt_allele_prob.txt",
        combined_aap,
        delimiter="\t",
        header="\t".join(sorted_mf),
        comments="",
    )


def write_pheno_penetrance(pedigree):
    """Writes out the phenotype penetrance for each individual in the pedigree.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :return: None. Writes to a file specified in the InputOutput.args.
    """
    args = InputOutput.args
    np.savetxt(
        args.out_file + ".pheno_penetrance.txt",
        pedigree.phenoPenetrance,
    )


def write_genotypes(pedigree, geno_prob_func, is_x_chr):
    """Writes out the genotypes for each individual. Format depends on user input and arguments.
    The output can include:
    - Dosages
    - Phased genotype probabilities
    - Genotype probabilities
    - Called genotypes (based on user inputted threshold or default 1/3)
    - Called haplotypes (based on user inputted threshold or default 1/2)

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param is_x_chr: flag whether inputted genotypes are X chromosome
    :type is_x_chr: bool
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    x_chr_male_dosage_weights = np.array([0, 0, 0, 1])
    autosome_dosage_weights = np.array([0, 1, 1, 2])

    geno_threshold_list = []
    if args.geno:
        if args.geno_threshold:
            for thresh in args.geno_threshold:
                if thresh < 1 / 3:
                    geno_threshold_list.append(1 / 3)
                else:
                    geno_threshold_list.append(thresh)
        else:
            geno_threshold_list.append(1 / 3)

    hap_threshold_list = []
    if args.hap:
        if args.hap_threshold:
            for thresh in args.hap_threshold:
                if thresh < 1 / 2:
                    hap_threshold_list.append(1 / 2)
                else:
                    hap_threshold_list.append(thresh)
        else:
            hap_threshold_list.append(1 / 2)

    has_output = (
        (not args.no_dosage)
        or args.phased_geno_prob
        or args.geno_prob
        or len(geno_threshold_list) != 0
        or len(hap_threshold_list) != 0
    )
    if not has_output:
        return

    with ExitStack() as stack:
        dosage_file = None
        if not args.no_dosage:
            dosage_file = stack.enter_context(open(args.out_file + ".dosage.txt", "w+"))

        phased_geno_prob_file = None
        if args.phased_geno_prob:
            phased_geno_prob_file = stack.enter_context(
                open(args.out_file + ".phased_geno_prob.txt", "w+")
            )

        geno_prob_file = None
        if args.geno_prob:
            geno_prob_file = stack.enter_context(
                open(args.out_file + ".geno_prob.txt", "w+")
            )

        geno_files = []
        for threshold in geno_threshold_list:
            output_file = args.out_file + ".geno_" + str(round(threshold, 3)) + ".txt"
            print(
                f"Writing called genotypes with threshold {threshold} to {output_file}"
            )
            geno_files.append((threshold, stack.enter_context(open(output_file, "w+"))))

        hap_files = []
        for threshold in hap_threshold_list:
            output_file = args.out_file + ".hap_" + str(round(threshold, 3)) + ".txt"
            print(
                f"Writing called haplotypes with threshold {threshold} to {output_file}"
            )
            hap_files.append((threshold, stack.enter_context(open(output_file, "w+"))))

        if is_x_chr:
            for idx, ind in pedigree.writeOrder():
                matrix = geno_prob_func(ind.idn, ind.sex)

                if dosage_file is not None:
                    write_dosage_from_matrix(
                        dosage_file,
                        ind,
                        matrix,
                        True,
                        formatter,
                        x_chr_male_dosage_weights,
                        autosome_dosage_weights,
                    )

                if phased_geno_prob_file is not None:
                    write_phased_geno_probs_from_matrix(
                        phased_geno_prob_file, ind, matrix, formatter
                    )

                if geno_prob_file is not None:
                    write_geno_probs_from_matrix(geno_prob_file, ind, matrix, formatter)

                if geno_files:
                    matrix_collapsed_hets = get_collapsed_genotypes(
                        matrix, True, ind.sex
                    )
                    for threshold, output_handle in geno_files:
                        write_called_genotypes_from_collapsed(
                            output_handle, ind, matrix_collapsed_hets, threshold
                        )

                if hap_files:
                    for threshold, output_handle in hap_files:
                        write_called_phase_from_matrix(
                            output_handle, ind, matrix, True, threshold
                        )
        else:
            for idx, ind in pedigree.writeOrder():
                matrix = geno_prob_func(ind.idn, ind.sex)

                if dosage_file is not None:
                    write_autosome_dosage_from_matrix(
                        dosage_file,
                        ind,
                        matrix,
                        formatter,
                        autosome_dosage_weights,
                    )

                if phased_geno_prob_file is not None:
                    write_phased_geno_probs_from_matrix(
                        phased_geno_prob_file, ind, matrix, formatter
                    )

                if geno_prob_file is not None:
                    write_geno_probs_from_matrix(geno_prob_file, ind, matrix, formatter)

                if geno_files:
                    matrix_collapsed_hets = get_autosome_collapsed_genotypes(matrix)
                    for threshold, output_handle in geno_files:
                        write_called_genotypes_from_collapsed(
                            output_handle, ind, matrix_collapsed_hets, threshold
                        )

                if hap_files:
                    for threshold, output_handle in hap_files:
                        write_autosome_called_phase_from_matrix(
                            output_handle, ind, matrix, threshold
                        )


def write_genotypes_separate_passes(pedigree, geno_prob_func, is_x_chr):
    """Writes genotype outputs with one full pedigree pass per output file."""

    args = InputOutput.args
    if not args.no_dosage:
        write_dosages(pedigree, geno_prob_func, is_x_chr, args.out_file + ".dosage.txt")
    if args.phased_geno_prob:
        write_phased_geno_probs(
            pedigree, geno_prob_func, args.out_file + ".phased_geno_prob.txt"
        )
    if args.geno_prob:
        write_geno_probs(pedigree, geno_prob_func, args.out_file + ".geno_prob.txt")
    if args.geno:
        geno_threshold_list = []
        if args.geno_threshold:
            for thresh in args.geno_threshold:
                if thresh < 1 / 3:
                    geno_threshold_list.append(1 / 3)
                else:
                    geno_threshold_list.append(thresh)
        else:
            geno_threshold_list.append(1 / 3)

        for threshold in geno_threshold_list:
            print(
                f"Writing called genotypes with threshold {threshold} to {args.out_file + '.geno_' + str(round(threshold, 3)) + '.txt'}"
            )
            write_called_genotypes(
                pedigree,
                geno_prob_func,
                is_x_chr,
                args.out_file + ".geno_" + str(round(threshold, 3)) + ".txt",
                threshold,
            )

    if args.hap:
        hap_threshold_list = []
        if args.hap_threshold:
            for thresh in args.hap_threshold:
                if thresh < 1 / 2:
                    hap_threshold_list.append(1 / 2)
                else:
                    hap_threshold_list.append(thresh)
        else:
            hap_threshold_list.append(1 / 2)

        for threshold in hap_threshold_list:
            print(
                f"Writing called haplotypes with threshold {threshold} to {args.out_file + '.hap_' + str(round(threshold, 3)) + '.txt'}"
            )
            write_called_phase(
                pedigree,
                geno_prob_func,
                is_x_chr,
                args.out_file + ".hap_" + str(round(threshold, 3)) + ".txt",
                threshold,
            )


def write_phased_geno_probs_from_matrix(f, ind, matrix, formatter):
    """Writes one individual's phased genotype probabilities."""

    for i in range(matrix.shape[0]):
        matrix_row = matrix[i, :]
        f.write(ind.idx + " " + " ".join(map(formatter, matrix_row)) + "\n")


def write_geno_probs_from_matrix(f, ind, matrix, formatter):
    """Writes one individual's unphased genotype probabilities."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]
    for i in range(matrix.shape[0]):
        if i == 1:  # Add up probabilities for aA and Aa
            f.write(ind.idx + " " + " ".join(map(formatter, matrix1 + matrix2)) + "\n")
        elif i != 2:  # Print probabilities for aa and AA
            matrix_row = matrix0
            if i == 3:
                matrix_row = matrix3
            f.write(ind.idx + " " + " ".join(map(formatter, matrix_row)) + "\n")


def write_dosage_from_matrix(
    f,
    ind,
    matrix,
    is_x_chr,
    formatter,
    x_chr_male_dosage_weights,
    autosome_dosage_weights,
):
    """Writes one individual's allele dosage."""

    if is_x_chr and ind.sex == 0:
        tmp = x_chr_male_dosage_weights
    else:
        tmp = autosome_dosage_weights
    dosage = np.dot(tmp, matrix)
    f.write(ind.idx + " " + " ".join(map(formatter, dosage)) + "\n")


def write_autosome_dosage_from_matrix(f, ind, matrix, formatter, dosage_weights):
    """Writes one individual's autosomal allele dosage."""

    dosage = np.dot(dosage_weights, matrix)
    f.write(ind.idx + " " + " ".join(map(formatter, dosage)) + "\n")


def get_collapsed_genotypes(matrix, is_x_chr, sex):
    """Collapse phased genotype probabilities into called-genotype states."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]
    if is_x_chr and sex == 0:
        matrix_collapsed_hets = np.empty((2, matrix.shape[1]), dtype=np.float32)
        matrix_collapsed_hets[0, :] = matrix0 + matrix2
        matrix_collapsed_hets[1, :] = matrix1 + matrix3
    else:
        matrix_collapsed_hets = np.empty((3, matrix.shape[1]), dtype=np.float32)
        matrix_collapsed_hets[0, :] = matrix0
        matrix_collapsed_hets[1, :] = matrix1 + matrix2
        matrix_collapsed_hets[2, :] = matrix3
    return matrix_collapsed_hets


def get_autosome_collapsed_genotypes(matrix):
    """Collapse autosomal phased genotype probabilities into called-genotype states."""

    matrix_collapsed_hets = np.empty((3, matrix.shape[1]), dtype=np.float32)
    matrix_collapsed_hets[0, :] = matrix[0, :]
    matrix_collapsed_hets[1, :] = matrix[1, :] + matrix[2, :]
    matrix_collapsed_hets[2, :] = matrix[3, :]
    return matrix_collapsed_hets


def write_called_genotypes_from_collapsed(f, ind, matrix_collapsed_hets, thresh):
    """Writes one individual's called genotypes from collapsed probabilities."""

    called_genotypes = np.argmax(matrix_collapsed_hets, axis=0)
    set_missing(called_genotypes, matrix_collapsed_hets, thresh)
    f.write(ind.idx + " " + " ".join(map(str, called_genotypes)) + "\n")


def write_called_phase_from_matrix(f, ind, matrix, is_x_chr, thresh):
    """Writes one individual's called haplotypes."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]

    if is_x_chr and ind.sex == 0:
        paternal_haplotype = np.full(matrix.shape[1], 9, dtype=np.int8)
    else:
        paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
        paternal_probs[0, :] = matrix0 + matrix1
        paternal_probs[1, :] = matrix2 + matrix3
        paternal_haplotype = np.argmax(paternal_probs, axis=0)
        set_missing(paternal_haplotype, paternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

    maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    maternal_probs[0, :] = matrix0 + matrix2
    maternal_probs[1, :] = matrix1 + matrix3
    maternal_haplotype = np.argmax(maternal_probs, axis=0)
    set_missing(maternal_haplotype, maternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


def write_autosome_called_phase_from_matrix(f, ind, matrix, thresh):
    """Writes one individual's autosomal called haplotypes."""

    matrix0 = matrix[0, :]
    matrix1 = matrix[1, :]
    matrix2 = matrix[2, :]
    matrix3 = matrix[3, :]

    paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    paternal_probs[0, :] = matrix0 + matrix1
    paternal_probs[1, :] = matrix2 + matrix3
    paternal_haplotype = np.argmax(paternal_probs, axis=0)
    set_missing(paternal_haplotype, paternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

    maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
    maternal_probs[0, :] = matrix0 + matrix2
    maternal_probs[1, :] = matrix1 + matrix3
    maternal_haplotype = np.argmax(maternal_probs, axis=0)
    set_missing(maternal_haplotype, maternal_probs, thresh)
    f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


def write_phased_geno_probs(pedigree, geno_prob_func, output_file):
    """Writes the phased genotype probabilities to a file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param output_file: name of output file to write to
    :type output_file: str
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(output_file, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            for i in range(matrix.shape[0]):
                matrix_row = matrix[i, :]
                f.write(ind.idx + " " + " ".join(map(formatter, matrix_row)) + "\n")


def write_geno_probs(pedigree, geno_prob_func, output_file):
    """Writes out the non phased genotype probabilities to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param output_file: name of output file to write to
    :type output_file: str
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(output_file, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]
            for i in range(matrix.shape[0]):
                if i == 1:  # Add up probabilities for aA and Aa
                    f.write(
                        ind.idx
                        + " "
                        + " ".join(
                            map(
                                formatter,
                                matrix1 + matrix2,
                            )
                        )
                        + "\n"
                    )
                elif i != 2:  # Print probabilities for aa and AA
                    matrix_row = matrix0
                    if i == 3:
                        matrix_row = matrix3
                    f.write(ind.idx + " " + " ".join(map(formatter, matrix_row)) + "\n")


def write_pheno_probs(pedigree, pheno_prob_func):
    """Writes out the phenotype probabilities to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(args.out_file + ".pheno_prob.txt", "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = pheno_prob_func(ind.idn, pedigree.phenoPenetrance)
            f.write("\n")
            for i in range(matrix.shape[0]):
                matrix_row = matrix[i, :]
                f.write(ind.idx + " " + " ".join(map(formatter, matrix_row)) + "\n")


def write_dosages(pedigree, geno_prob_func, is_x_chr, output_file):
    """Writes out the allele dosages to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param is_x_chr: flag whether inputted genotypes are X chromosome
    :type is_x_chr: bool
    :param output_file: name of output file to write to
    :type output_file: str
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    x_chr_male_dosage_weights = np.array([0, 0, 0, 1])
    autosome_dosage_weights = np.array([0, 1, 1, 2])
    with open(output_file, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            if is_x_chr and ind.sex == 0:
                tmp = x_chr_male_dosage_weights
            else:
                tmp = autosome_dosage_weights
            matrix = np.dot(tmp, geno_prob_func(ind.idn, ind.sex))
            f.write(ind.idx + " " + " ".join(map(formatter, matrix)) + "\n")


def write_called_genotypes(pedigree, geno_prob_func, is_x_chr, output_file, thresh):
    """Writes out the called genotypes to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param is_x_chr: flag whether inputted genotypes are X chromosome
    :type is_x_chr: bool
    :param output_file: name of output file to write to
    :type output_file: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(output_file, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]
            if is_x_chr and ind.sex == 0:
                matrix_collapsed_hets = np.empty((2, matrix.shape[1]), dtype=np.float32)
                matrix_collapsed_hets[0, :] = matrix0 + matrix2
                matrix_collapsed_hets[1, :] = matrix1 + matrix3
            else:
                matrix_collapsed_hets = np.empty((3, matrix.shape[1]), dtype=np.float32)
                matrix_collapsed_hets[0, :] = matrix0
                matrix_collapsed_hets[1, :] = matrix1 + matrix2
                matrix_collapsed_hets[2, :] = matrix3

            called_genotypes = np.argmax(matrix_collapsed_hets, axis=0)
            set_missing(called_genotypes, matrix_collapsed_hets, thresh)
            f.write(ind.idx + " " + " ".join(map(str, called_genotypes)) + "\n")


def write_called_phase(pedigree, geno_prob_func, is_x_chr, output_file, thresh):
    """Writes out the called haplotypes to file.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param output_file: name of output file to write to
    :type output_file: str
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(output_file, "w+") as f:
        for idx, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            matrix0 = matrix[0, :]
            matrix1 = matrix[1, :]
            matrix2 = matrix[2, :]
            matrix3 = matrix[3, :]

            # Paternal
            if is_x_chr and ind.sex == 0:
                paternal_haplotype = np.full(matrix.shape[1], 9, dtype=np.int8)
            else:
                paternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
                paternal_probs[0, :] = matrix0 + matrix1
                paternal_probs[1, :] = matrix2 + matrix3
                paternal_haplotype = np.argmax(paternal_probs, axis=0)
                set_missing(paternal_haplotype, paternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, paternal_haplotype)) + "\n")

            # Maternal
            maternal_probs = np.empty((2, matrix.shape[1]), dtype=np.float32)
            maternal_probs[0, :] = matrix0 + matrix2
            maternal_probs[1, :] = matrix1 + matrix3
            maternal_haplotype = np.argmax(maternal_probs, axis=0)
            set_missing(maternal_haplotype, maternal_probs, thresh)
            f.write(ind.idx + " " + " ".join(map(str, maternal_haplotype)) + "\n")


@jit(nopython=True)
def set_missing(called_genotypes, matrix, thresh):
    """Sets the called genotypes to missing if the probability is below the threshold.

    :param called_genotypes: array of called genotypes
    :type called_genotypes: 2D numpy array of float32 with shape 3 x n_loci
    :param matrix: genotype probability matrix
    :type matrix: 2D numpy array of float32 with shape 3 x n_loci
    :param thresh: threshold for calling genotypes, defaults to 1/3
    :type thresh: float
    """
    n_loci = len(called_genotypes)
    for i in range(n_loci):
        if matrix[called_genotypes[i], i] <= thresh:
            called_genotypes[i] = 9
