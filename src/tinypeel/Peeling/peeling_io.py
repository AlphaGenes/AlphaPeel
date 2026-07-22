"""Input/output functions."""

from contextlib import ExitStack

import numpy as np
from numba import jit
from ..tinyhouse import InputOutput


def write_out_parameters(peeling_info):
    """Write estimated error rates and recombination probabilities.

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
        # Placeholder until recombination-probability output is implemented.
        np.savetxt(
            args.out_file + ".rec_prob.txt",
            np.empty((1, 1)),
        )


def write_out_alt_allele_prob(pedigree):
    """Write alternative allele probabilities for each locus and metafounder.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :return: None. Writes to a file specified in the InputOutput.args.
    """
    args = InputOutput.args

    def sort_key(mf_key):
        part = mf_key.split("_")[1]
        return (0, int(part)) if part.isdigit() else (1, part)

    # Sort MF_1, MF_2, ..., MF_11 numerically when possible.
    sorted_aap = dict(sorted(pedigree.AAP.items(), key=lambda item: sort_key(item[0])))
    sorted_mf = list(sorted_aap.keys())
    combined_aap = np.hstack(
        [sorted_aap[key].reshape(pedigree.nLoci, -1) for key in sorted_mf]
    )
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
    """Write requested genotype outputs for each individual.

    Depending on CLI options, outputs can include dosages, phased genotype
    probabilities, genotype probabilities, called genotypes, and called
    haplotypes.

    :param pedigree: pedigree information container
    :type pedigree: class:`tinyhouse.Pedigree.Pedigree()`
    :param geno_prob_func: function to get genotype probabilities for an individual
    :type geno_prob_func: function
    :param is_x_chr: flag whether inputted genotypes are X chromosome
    :type is_x_chr: bool
    :return: None. Writes to files specified in the InputOutput.args.
    """
    args = InputOutput.args
    geno_thresholds = get_output_thresholds(args.geno, args.geno_threshold, 1 / 3)
    hap_thresholds = get_output_thresholds(args.hap, args.hap_threshold, 1 / 2)
    output_settings = {
        "formatter": f"{{:.{args.out_digits}f}}".format,
        "dosage_weights": (np.array([0, 0, 0, 1]), np.array([0, 1, 1, 2])),
        "is_x_chr": is_x_chr,
        "geno_thresholds": geno_thresholds,
        "hap_thresholds": hap_thresholds,
    }

    if not (
        (not args.no_dosage)
        or args.phased_geno_prob
        or args.geno_prob
        or len(geno_thresholds) != 0
        or len(hap_thresholds) != 0
    ):
        return

    write_genotypes_single_pass(pedigree, geno_prob_func, output_settings)


def get_output_thresholds(enabled, thresholds, minimum):
    """Return validated output thresholds for called genotypes or haplotypes."""

    if not enabled:
        return []
    if thresholds:
        return [max(thresh, minimum) for thresh in thresholds]
    return [minimum]


def write_genotypes_single_pass(pedigree, geno_prob_func, output_settings):
    """Write all requested genotype outputs in one pass through the pedigree."""

    args = InputOutput.args
    with ExitStack() as stack:
        output_files = open_genotype_output_files(
            stack, args, output_settings["geno_thresholds"]
        )
        output_files["haps"] = open_threshold_files(
            stack, args.out_file, "hap", output_settings["hap_thresholds"], "haplotypes"
        )

        for _, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            write_individual_genotype_outputs(
                ind, matrix, output_settings, output_files
            )


def open_genotype_output_files(stack, args, geno_thresholds):
    """Open requested genotype output files and return their handles."""

    output_files = {
        "dosage": None,
        "phased_geno_prob": None,
        "geno_prob": None,
        "genotypes": open_threshold_files(
            stack, args.out_file, "geno", geno_thresholds, "genotypes"
        ),
    }
    if not args.no_dosage:
        output_files["dosage"] = enter_output_file(stack, args.out_file + ".dosage.txt")
    if args.phased_geno_prob:
        output_files["phased_geno_prob"] = enter_output_file(
            stack, args.out_file + ".phased_geno_prob.txt"
        )
    if args.geno_prob:
        output_files["geno_prob"] = enter_output_file(
            stack, args.out_file + ".geno_prob.txt"
        )
    return output_files


def open_threshold_files(stack, out_file, file_suffix, thresholds, output_label):
    """Open requested thresholded output files."""

    output_files = []
    for threshold in thresholds:
        output_file = f"{out_file}.{file_suffix}_{round(threshold, 3)}.txt"
        print(
            f"Writing called {output_label} with threshold {threshold} to {output_file}"
        )
        output_files.append((threshold, enter_output_file(stack, output_file)))
    return output_files


def enter_output_file(stack, output_file):
    """Open an output file under the shared output ExitStack."""

    # pylint: disable=consider-using-with
    return stack.enter_context(open(output_file, "w+", encoding="utf-8"))


def write_individual_genotype_outputs(ind, matrix, output_settings, output_files):
    """Write all requested genotype outputs for one individual."""

    if output_files["dosage"] is not None:
        write_individual_dosage(output_files["dosage"], ind, matrix, output_settings)
    if output_files["phased_geno_prob"] is not None:
        write_phased_geno_probs_from_matrix(
            output_files["phased_geno_prob"],
            ind,
            matrix,
            output_settings["formatter"],
        )
    if output_files["geno_prob"] is not None:
        write_geno_probs_from_matrix(
            output_files["geno_prob"], ind, matrix, output_settings["formatter"]
        )
    if output_files["genotypes"]:
        write_individual_called_genotypes(ind, matrix, output_settings, output_files)
    if output_files["haps"]:
        write_individual_called_haplotypes(
            ind, matrix, output_settings["is_x_chr"], output_files["haps"]
        )


def write_individual_dosage(output_file, ind, matrix, output_settings):
    """Write one individual's dosage to the appropriate chromosome output."""

    x_chr_male_weights, autosome_weights = output_settings["dosage_weights"]
    if output_settings["is_x_chr"]:
        write_dosage_from_matrix(
            output_file,
            ind,
            matrix,
            output_settings["formatter"],
            (x_chr_male_weights, autosome_weights),
        )
    else:
        write_autosome_dosage_from_matrix(
            output_file,
            ind,
            matrix,
            output_settings["formatter"],
            autosome_weights,
        )


def write_individual_called_genotypes(ind, matrix, output_settings, output_files):
    """Write one individual's called genotype outputs."""

    if output_settings["is_x_chr"]:
        matrix_collapsed_hets = get_collapsed_genotypes(matrix, True, ind.sex)
    else:
        matrix_collapsed_hets = get_autosome_collapsed_genotypes(matrix)
    for threshold, output_handle in output_files["genotypes"]:
        write_called_genotypes_from_collapsed(
            output_handle, ind, matrix_collapsed_hets, threshold
        )


def write_individual_called_haplotypes(ind, matrix, is_x_chr, hap_files):
    """Write one individual's called haplotype outputs."""

    for threshold, output_handle in hap_files:
        if is_x_chr:
            write_called_phase_from_matrix(output_handle, ind, matrix, True, threshold)
        else:
            write_autosome_called_phase_from_matrix(
                output_handle, ind, matrix, threshold
            )


def write_genotypes_separate_passes(pedigree, geno_prob_func, is_x_chr):
    """Writes genotype outputs with one full pedigree pass per output file."""

    args = InputOutput.args
    output_context = {
        "pedigree": pedigree,
        "geno_prob_func": geno_prob_func,
        "is_x_chr": is_x_chr,
        "out_file": args.out_file,
    }
    if not args.no_dosage:
        write_dosages(pedigree, geno_prob_func, is_x_chr, args.out_file + ".dosage.txt")
    if args.phased_geno_prob:
        write_phased_geno_probs(
            pedigree, geno_prob_func, args.out_file + ".phased_geno_prob.txt"
        )
    if args.geno_prob:
        write_geno_probs(pedigree, geno_prob_func, args.out_file + ".geno_prob.txt")
    write_separate_threshold_outputs(
        output_context,
        get_separate_output_settings(
            args.geno, args.geno_threshold, 1 / 3, "geno", "genotypes"
        ),
        write_called_genotypes,
    )
    write_separate_threshold_outputs(
        output_context,
        get_separate_output_settings(
            args.hap, args.hap_threshold, 1 / 2, "hap", "haplotypes"
        ),
        write_called_phase,
    )


def get_separate_output_settings(enabled, thresholds, minimum, suffix, label):
    """Return settings for thresholded separate-pass outputs."""

    return {
        "thresholds": get_output_thresholds(enabled, thresholds, minimum),
        "suffix": suffix,
        "label": label,
    }


def write_separate_threshold_outputs(output_context, output_settings, writer):
    """Write separate-pass thresholded genotype or haplotype outputs."""

    for threshold in output_settings["thresholds"]:
        output_file = (
            f"{output_context['out_file']}.{output_settings['suffix']}_"
            f"{round(threshold, 3)}.txt"
        )
        print(
            f"Writing called {output_settings['label']} "
            f"with threshold {threshold} to {output_file}"
        )
        writer(
            output_context["pedigree"],
            output_context["geno_prob_func"],
            output_context["is_x_chr"],
            output_file,
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
    formatter,
    dosage_weights,
):
    """Writes one individual's allele dosage."""

    x_chr_male_dosage_weights, autosome_dosage_weights = dosage_weights
    if ind.sex == 0:
        weights = x_chr_male_dosage_weights
    else:
        weights = autosome_dosage_weights
    dosage = np.dot(weights, matrix)
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
    with open(output_file, "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
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
    with open(output_file, "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
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
    :param pheno_prob_func: function to get phenotype probabilities for an individual
    :type pheno_prob_func: function
    :return: None. Writes to the specified output file.
    """
    args = InputOutput.args
    formatter = f"{{:.{args.out_digits}f}}".format
    with open(args.out_file + ".pheno_prob.txt", "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
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
    with open(output_file, "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
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
    with open(output_file, "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
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
    :param is_x_chr: flag whether inputted genotypes are X chromosome
    :type is_x_chr: bool
    :param output_file: name of output file to write to
    :type output_file: str
    :param thresh: threshold for calling haplotypes, defaults to 1/2
    :type thresh: float
    :return: None. Writes to the specified output file.
    """
    with open(output_file, "w+", encoding="utf-8") as f:
        for _, ind in pedigree.writeOrder():
            matrix = geno_prob_func(ind.idn, ind.sex)
            write_called_phase_from_matrix(f, ind, matrix, is_x_chr, thresh)


@jit(nopython=True)
def set_missing(called_types, matrix, thresh):
    """Sets the called genotypes or haplotypes to missing if the probability is below the threshold.

    :param called_types: array of called genotypes or haplotypes
    :type called_types: 1D numpy array of integer genotype or haplotype calls
    :param matrix: genotype probability matrix
    :type matrix: 2D numpy array of probabilities with states x n_loci
    :param thresh: threshold for calling genotypes or haplotypes
    :type thresh: float
    """
    n_loci = len(called_types)
    for i in range(n_loci):
        if matrix[called_types[i], i] <= thresh:
            called_types[i] = 9
