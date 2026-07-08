import numpy as np
import os
import subprocess
import argparse
import warnings
import shutil
import src.tinypeel.tinypeel as tinypeel


def get_params():
    param_file = os.path.join("tests", "accuracy_tests", "simulation_parameters.txt")
    with open(param_file, "r") as file:
        sim_params = [line.strip().split() for line in file]

    params = {}
    for param_name, param_value in sim_params:
        params[param_name] = float(param_value)

    return params


def sim_path():
    """Return the simulation directory used by the accuracy tests.

    :return: Relative path to the simulation fixture directory.
    :rtype: str
    """

    return os.path.join("tests", "accuracy_tests", "sim_for_alphapeel_accu_test")


def get_accuracy_benchmark_output_root(run_name="test_accu"):
    """Return the root directory for direct-call accuracy benchmark outputs."""

    return os.path.join("tests", "accuracy_tests", f"outputs_{run_name}")


def get_accuracy_benchmark_report_root(run_name="test_accu"):
    """Return the root directory for direct-call accuracy benchmark reports."""

    return os.path.join("tests", "accuracy_tests", f"reports_{run_name}")


def prepare_directory(path):
    """Create an empty directory at ``path``."""

    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path, exist_ok=True)


def ensure_directory(path):
    """Create ``path`` if it does not exist."""

    os.makedirs(path, exist_ok=True)


def build_accuracy_output_path(name, run_name="test_accu"):
    """Build the per-case direct-call benchmark output directory path."""

    return os.path.join(get_accuracy_benchmark_output_root(run_name), name)


def build_accuracy_report_path(run_name="test_accu"):
    """Build the direct-call benchmark report file path."""

    return os.path.join(get_accuracy_benchmark_report_root(run_name), "accu_report.txt")


def run_command(command):
    result = subprocess.run(
        command,
        shell=True,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "AlphaPeel run failed\n"
            f"Command: {command}\n"
            f"Exit code: {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    return result


def get_accuracy_benchmark_cases():
    """Return the full accuracy benchmark case matrix."""

    return [
        ("single", None, None, None, None, None, None, None, None),
        (
            "single",
            "est_start_alt_allele_prob",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ),
        ("single", None, None, None, None, None, "est_alt_allele_prob", None, None),
        ("multi", None, None, None, None, None, None, None, None),
        (
            "multi",
            "est_start_alt_allele_prob",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ),
        ("multi", None, None, None, None, None, "est_alt_allele_prob", None, None),
        (
            "multi",
            "est_start_alt_allele_prob",
            "est_geno_error_prob",
            "est_seq_error_prob",
            None,
            None,
            None,
            None,
            None,
        ),
        ("multi", None, None, None, "seq_file", None, None, None, None),
        (
            "multi",
            "est_start_alt_allele_prob",
            None,
            None,
            "seq_file",
            None,
            None,
            None,
            None,
        ),
        (
            "multi",
            "est_start_alt_allele_prob",
            "est_geno_error_prob",
            "est_seq_error_prob",
            "seq_file",
            None,
            None,
            None,
            None,
        ),
        ("hybrid", None, None, None, None, None, None, None, None),
        ("hybrid", None, None, None, "seq_file", None, None, None, None),
        (
            "single",
            None,
            None,
            None,
            None,
            "alt_allele_prob_file",
            None,
            "metafounder",
            None,
        ),
        (
            "single",
            "est_start_alt_allele_prob",
            None,
            None,
            None,
            None,
            None,
            "metafounder",
            None,
        ),
        (
            "single",
            None,
            None,
            None,
            None,
            None,
            "est_alt_allele_prob",
            "metafounder",
            None,
        ),
        (
            "single",
            "est_start_alt_allele_prob",
            None,
            None,
            None,
            None,
            "est_alt_allele_prob",
            "metafounder",
            None,
        ),
        (
            "single",
            None,
            None,
            None,
            None,
            "alt_allele_prob_file",
            "est_alt_allele_prob",
            "metafounder",
            None,
        ),
        ("single", None, None, None, None, None, None, None, "x_chr"),
        ("multi", None, None, None, None, None, None, None, "x_chr"),
        ("hybrid", None, None, None, None, None, None, None, "x_chr"),
    ]


def build_accuracy_case_name(
    method,
    est_start_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    est_alt_allele_prob,
    metafounder,
    x_chr,
):
    """Build the benchmark case name from the accuracy case parameters."""

    return "_".join(
        [
            param
            for param in filter(
                lambda param: True if param else False,
                [
                    method,
                    est_start_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    est_alt_allele_prob,
                    metafounder,
                    x_chr,
                ],
            )
        ]
    )


def generate_accuracy_argv(
    sim_path,
    method,
    est_start_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    est_alt_allele_prob,
    metafounder,
    x_chr,
    output_path,
):
    """Generate an ``argv`` list for a direct ``tinypeel.main`` accuracy run."""

    argv = []
    input_files = ["ped_file"]
    arguments = {
        "method": method,
        "n_cycle": "5",
        "n_thread": "6",
        "geno_threshold": ".1",
        "hap_threshold": ".1",
        "geno": None,
        "hap": None,
        "seg_prob": None,
        "geno_prob": None,
        "phased_geno_prob": None,
    }

    if est_start_alt_allele_prob:
        arguments["est_start_alt_allele_prob"] = None
    if est_geno_error_prob and est_seq_error_prob:
        arguments["est_geno_error_prob"] = None
        arguments["est_seq_error_prob"] = None
    if est_alt_allele_prob:
        arguments["est_alt_allele_prob"] = None
    if seq_file:
        input_files.append("seq_file")
    else:
        input_files.append("geno_file")
    if alt_allele_prob_file:
        input_files.append("alt_allele_prob_file")
    if method == "hybrid":
        input_files.append("map_file")
        input_files.append("seg_map_file")
        input_files.append("seg_file")

    if metafounder:
        for file_name in input_files:
            argv.extend(
                [
                    f"-{file_name}",
                    os.path.join(sim_path, f"metafounder_{file_name}.txt"),
                ]
            )
    elif x_chr:
        for file_name in input_files:
            argv.extend(
                [
                    f"-{file_name}",
                    os.path.join(sim_path, f"X_chr_{file_name}.txt"),
                ]
            )
        argv.append("-x_chr")
    else:
        for file_name in input_files:
            argv.extend([f"-{file_name}", os.path.join(sim_path, f"{file_name}.txt")])

    for key, value in arguments.items():
        argv.append(f"-{key}")
        if value is not None:
            argv.append(value)

    argv.extend(["-out_file", f"{output_path}{os.sep}"])

    return argv


def generate_command(
    sim_path,
    method,
    output_path,
):
    """Generate the shell command used to run AlphaPeel for a test case.

    :param sim_path: Directory containing simulation input files.
    :type sim_path: str
    :param method: AlphaPeel method to run, such as ``single``, ``multi`` or
        ``hybrid``.
    :type method: str
    :param output_path: Output directory for AlphaPeel results.
    :type output_path: str
    :return: The command string to pass to ``os.system``.
    :rtype: str
    """

    command = "AlphaPeel "
    input_file = ["ped_file"]
    arguments = {
        "method": method,
        "n_cycle": "5",
        "n_thread": "6",
        "geno_threshold": ".1",
        "hap_threshold": ".1",
        "geno": None,
        "hap": None,
        "seg_prob": None,
        "geno_prob": None,
        "phased_geno_prob": None,
    }

    input_file.append("geno_file")

    for file in input_file:
        command += f"-{file} {os.path.join(sim_path, f'{file}.txt')} "

    for key, value in arguments.items():
        if value is not None:
            command += f"-{key} {value} "
        else:
            command += f"-{key} "

    command += f"-out_file {output_path}{os.sep}"

    return command


def run_tinypeel_direct(argv):
    """Run AlphaPeel by calling ``src.tinypeel.tinypeel.main`` directly."""

    tinypeel.main(argv=argv)


def assess_test_accuracy(
    sim_path,
    get_params,
    output_path,
    method,
    file_out,
):
    file_to_check = [
        "dosage",
        "geno_0.333",
        "hap_0.5",
        "geno_prob",
        "phased_geno_prob",
    ]
    if method == "multi":
        file_to_check.append("seg_prob")

    nGen = int(get_params["nGen"])
    nIndPerGen = int(get_params["nInd"] / nGen)
    nLociAll = int(get_params["nLociAll"])

    for file_name in file_to_check:
        if file_name in ["dosage", "geno_0.333"]:
            n_row_per_ind = 1
        elif file_name in ["phased_geno_prob", "seg_prob"]:
            n_row_per_ind = 4
        elif file_name == "geno_prob":
            n_row_per_ind = 3
        elif file_name == "hap_0.5":
            n_row_per_ind = 2

        file_path = os.path.join(output_path, f".{file_name}.txt")
        true_path = os.path.join(sim_path, f"true-{file_name}.txt")

        try:
            new_file = np.loadtxt(file_path, usecols=np.arange(1, nLociAll + 1))
        except ValueError:
            print(f"Error loading {file_path}")
            continue

        try:
            true_file = np.loadtxt(true_path, usecols=np.arange(1, nLociAll + 1))
        except ValueError:
            print(f"Error loading {true_path}")
            continue

        if file_name == "seg_prob":
            marker_corr = get_marker_corr(
                new_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
                true_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
            )
        else:
            marker_corr = get_marker_corr(new_file[:, 1:], true_file[:, 1:])

        file_out.write(f"{file_name},{method},marker_corr,{marker_corr}\n")

        if file_name == "seg_prob":
            ind_corr = get_ind_corr(
                new_file[:, 1:],
                true_file[:, 1:],
                nIndPerGen,
                n_row_per_ind,
                start_gen=2,
            )
        else:
            ind_corr = get_ind_corr(
                new_file[:, 1:],
                true_file[:, 1:],
                nIndPerGen,
                n_row_per_ind,
            )

        file_out.write(f"{file_name},{method},ind_corr,{ind_corr}\n")

        if file_name == "segregation":
            abs_diff = get_abs_diff(
                new_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
                true_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
                n_row_per_ind,
            )

        else:
            abs_diff = get_abs_diff(
                new_file[:, 1:],
                true_file[:, 1:],
                n_row_per_ind,
            )

        file_out.write(f"{file_name},{method},abs_diff,{abs_diff}\n")

        if file_name == "segregation":
            correct_rate = get_correct_rate(
                new_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
                true_file[2 * (nIndPerGen * n_row_per_ind) :, 1:],
            )
            file_out.write(f"{file_name},{method},correct_rate,{correct_rate}\n")


def assess_accuracy(
    sim_path,
    get_params,
    output_path,
    name,
    method,
    file_out,
    metafounder,
    x_chr,
):
    """Assess output accuracy against truth files."""

    if metafounder:
        file_to_check = [
            "dosage",
            "geno_prob",
            "phased_geno_prob",
        ]
    else:
        file_to_check = [
            "dosage",
            "geno_0.333",
            "hap_0.5",
            "geno_prob",
            "phased_geno_prob",
        ]
    if method == "multi":
        file_to_check.append("seg_prob")

    nGen = int(get_params["nGen"])
    nIndPerGen = int(get_params["nInd"] / nGen)
    nLociAll = int(get_params["nLociAll"])

    print(" ")
    print(f"Test: {name}")

    for file_name in file_to_check:
        if file_name in ["dosage", "geno_0.333"]:
            n_row_per_ind = 1
        elif file_name in ["phased_geno_prob", "seg_prob"]:
            n_row_per_ind = 4
        elif file_name == "geno_prob":
            n_row_per_ind = 3
        elif file_name == "hap_0.5":
            n_row_per_ind = 2

        file_path = os.path.join(output_path, f".{file_name}.txt")
        if metafounder:
            true_path = os.path.join(sim_path, f"true-metafounder_{file_name}.txt")
        elif x_chr:
            true_path = os.path.join(sim_path, f"true-X_chr_{file_name}.txt")
        else:
            true_path = os.path.join(sim_path, f"true-{file_name}.txt")

        new_file = np.loadtxt(file_path, usecols=np.arange(1, nLociAll + 1))
        true_file = np.loadtxt(true_path, usecols=np.arange(1, nLociAll + 1))

        marker_corr = [str(get_marker_corr(new_file[:, :], true_file[:, :]))]
        for gen in range(nGen):
            marker_corr.append(
                str(
                    get_marker_corr(
                        new_file[
                            gen
                            * (nIndPerGen * n_row_per_ind) : (gen + 1)
                            * (nIndPerGen * n_row_per_ind)
                        ],
                        true_file[
                            gen
                            * (nIndPerGen * n_row_per_ind) : (gen + 1)
                            * (nIndPerGen * n_row_per_ind)
                        ],
                    )
                )
            )

        file_out.write(f"{file_name},{method},marker_corr,{marker_corr}\n")

        ind_corr = [
            str(get_ind_corr(new_file[:, :], true_file[:, :], nIndPerGen, None))
        ]
        for gen in range(nGen):
            ind_corr.append(
                str(
                    get_ind_corr(
                        new_file[:, 1:],
                        true_file[:, 1:],
                        nIndPerGen,
                        n_row_per_ind,
                        gen,
                        gen + 1,
                    )
                )
            )

        file_out.write(f"{file_name},{method},ind_corr,{ind_corr}\n")


def _prepare_hybrid_seg_file(
    get_params,
    sim_path,
    run_name,
    est_start_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    est_alt_allele_prob,
    metafounder,
    x_chr,
):
    """Create the hybrid segregation fixture from the matching multi output."""

    multi_name = build_accuracy_case_name(
        "multi",
        est_start_alt_allele_prob,
        est_geno_error_prob,
        est_seq_error_prob,
        seq_file,
        alt_allele_prob_file,
        est_alt_allele_prob,
        metafounder,
        x_chr,
    )
    multi_path = build_accuracy_output_path(multi_name, run_name)

    nSegMap = int(get_params["nSegMap"])
    nLociAll = int(get_params["nLociAll"])

    subset = np.floor(np.linspace(1, nLociAll, num=nSegMap)).astype(dtype=int)
    subset = np.concatenate(([0], subset))
    seg_path = os.path.join(multi_path, ".seg_prob.txt")
    seg_file_path = os.path.join(sim_path, "seg_file.txt")

    seg = np.loadtxt(seg_path)
    np.savetxt(seg_file_path, seg[:, subset])


def run_accuracy_case(
    get_params,
    sim_path,
    method,
    est_start_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    est_alt_allele_prob,
    metafounder,
    x_chr,
    benchmark=None,
    run_name="test_accu",
):
    """Run one accuracy benchmark case through ``tinypeel.main`` directly."""

    name = build_accuracy_case_name(
        method,
        est_start_alt_allele_prob,
        est_geno_error_prob,
        est_seq_error_prob,
        seq_file,
        alt_allele_prob_file,
        est_alt_allele_prob,
        metafounder,
        x_chr,
    )
    output_path = build_accuracy_output_path(name, run_name)
    prepare_directory(output_path)

    if method == "hybrid":
        _prepare_hybrid_seg_file(
            get_params,
            sim_path,
            run_name,
            est_start_alt_allele_prob,
            est_geno_error_prob,
            est_seq_error_prob,
            seq_file,
            alt_allele_prob_file,
            est_alt_allele_prob,
            metafounder,
            x_chr,
        )

    if run_name == "test_accu":
        command = generate_command(sim_path, method, output_path)
        benchmark(run_command, command)

    else:
        argv = generate_accuracy_argv(
            sim_path,
            method,
            est_start_alt_allele_prob,
            est_geno_error_prob,
            est_seq_error_prob,
            seq_file,
            alt_allele_prob_file,
            est_alt_allele_prob,
            metafounder,
            x_chr,
            output_path,
        )

        run_tinypeel_direct(argv)

    report_path = build_accuracy_report_path(run_name)
    ensure_directory(os.path.dirname(report_path))
    with open(report_path, "a") as file_out:
        if run_name == "test_accu":
            assess_test_accuracy(
                sim_path,
                get_params,
                output_path,
                method,
                file_out,
            )
        elif run_name == "benchmark":
            assess_accuracy(
                sim_path,
                get_params,
                output_path,
                name,
                method,
                file_out,
                metafounder,
                x_chr,
            )


def run_full_accuracy_suite(run_name="benchmark"):
    """Run the full direct-call accuracy benchmark suite."""

    prepare_directory(get_accuracy_benchmark_output_root(run_name))
    prepare_directory(get_accuracy_benchmark_report_root(run_name))

    for case in get_accuracy_benchmark_cases():
        run_accuracy_case(
            get_params(),
            sim_path(),
            *case,
            benchmark=None,
            run_name=run_name,
        )


def get_marker_corr(output, real):
    """Compute the average marker-wise Pearson correlation coefficient between two arrays.

    :param output: Output matrix to compare.
    :type output: numpy.ndarray
    :param real: Reference matrix to compare against.
    :type real: numpy.ndarray
    :return: Rounded mean marker correlation coefficient.
    :rtype: float
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        accus = np.array(
            [np.corrcoef(real[:, i], output[:, i])[0, 1] for i in range(real.shape[1])]
        )
        return round(np.nanmean(accus), 4)


def get_ind_corr(output, real, nIndPerGen, n_row_per_ind, start_gen=None, end_gen=None):
    """Compute the average individual-wise Pearson correlation coefficient between two arrays.

    :param output: Output matrix to compare.
    :type output: numpy.ndarray
    :param real: Reference matrix to compare against.
    :type real: numpy.ndarray
    :param nIndPerGen: Number of individuals per generation.
    :type nIndPerGen: int
    :param n_row_per_ind: Number of rows per individual in the output.
    :type n_row_per_ind: int or None
    :param start_gen: Optional starting generation index used to slice the correlation vector.
    :type start_gen: int or None
    :param end_gen: Optional ending generation index used to slice the correlation vector.
    :type end_gen: int or None
    :return: Rounded mean individual accuracy.
    :rtype: float
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        accus = np.array(
            [np.corrcoef(real[i, :], output[i, :])[0, 1] for i in range(real.shape[0])]
        )
        if type(start_gen) == int:
            if type(end_gen) != int:
                accus = accus[start_gen * (nIndPerGen * n_row_per_ind) :]
            else:
                accus = accus[
                    start_gen
                    * (nIndPerGen * n_row_per_ind) : end_gen
                    * (nIndPerGen * n_row_per_ind)
                ]

        return round(np.nanmean(accus), 4)


def get_abs_diff(output, real, n_row_per_ind):
    """Sum of absolute difference divided by the sum of the number of loci being counted"""
    return n_row_per_ind * np.sum(np.abs(output - real)) / real.size


def get_correct_rate(output, real):
    """
    Summing up the probabilities of the true state from the output data
    divided by the number of loci being counted

    :param output: Output data
    :type output: ndarray
    :param real: Real simulated data
    :type real: ndarray
    """
    return np.sum(output[real == 1]) / (np.size(real) / 4)


def calc_switch_error_rate():
    """
    Run this function at root directory and input the path to the
    assessed haplotype file.
    The SER and PER_intra calculation follows the definition from:
    https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X
    """

    parser = argparse.ArgumentParser(
        prog="calc_switch_error_rate",
        description="Calculate switch error rate, phase error rate",
        epilog="The SER and PER_intra calculation follows the definition from: https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X",
    )

    params = get_params()

    nLociAll = int(params["nLociAll"])
    nInd = int(params["nInd"])

    parser.add_argument(
        "-true_path",
        type=str,
        required=False,
        default=os.path.join(
            "tests", "accuracy_tests", "sim_for_alphapeel_accu_test", "true-hap_0.5.txt"
        ),
        help="Enter the path of the true haplotype file, default is the simulation path",
    )
    parser.add_argument(
        "-called_path",
        type=str,
        required=True,
        help="Enter the path of the assessed haplotype file",
    )

    args = parser.parse_args()
    true_path = args.true_path
    called_path = args.called_path

    called_file = np.loadtxt(called_path, usecols=np.arange(1, nLociAll + 1))
    true_file = np.loadtxt(true_path, usecols=np.arange(1, nLociAll + 1))

    switch_error_count = 0
    phase_error_count = 0
    uncalled_count = 0
    wrong_homo_count = 0
    true_hetero_count = 0
    homo_count = 0
    hetero_count = 0

    for ind in range(nInd):
        hap_p_new = called_file[ind * 2]
        hap_m_new = called_file[ind * 2 + 1]
        hap_p_true = true_file[ind * 2]
        hap_m_true = true_file[ind * 2 + 1]

        switched = False
        for loci in range(nLociAll):
            if hap_p_true[loci] == hap_m_true[loci]:
                homo_count += 1
            else:
                hetero_count += 1
            if hap_p_new[loci] == 9 or hap_m_new[loci] == 9:
                uncalled_count += 1
                continue
            if hap_p_true[loci] + hap_m_true[loci] == 1:
                if hap_p_new[loci] == hap_m_new[loci]:
                    wrong_homo_count += 1
                    continue
                true_hetero_count += 1
                if hap_p_new[loci] != hap_p_true[loci]:
                    phase_error_count += 1
                if (hap_p_new[loci] != hap_p_true[loci] and switched is False) or (
                    hap_p_new[loci] == hap_p_true[loci] and switched is True
                ):
                    switched = not switched
                    switch_error_count += 1

    print(f"Switch error rate: {switch_error_count / (nInd * (nLociAll - 1))}")
    print(f"Phase error (intra) rate: {phase_error_count / (nInd * nLociAll)}")
    print(f"Uncalled rate: {uncalled_count / (nInd * nLociAll)}")
    print(
        f"Proportion of genotypes wrongly called as homozygote: {wrong_homo_count / (nInd * nLociAll)}"
    )
    print(
        f"Proportion of genotypes correctly called as heterozygote: {true_hetero_count / (nInd * nLociAll)}"
    )
    print(f"Homozygote count in true genotype: {homo_count}")
    print(f"Heterozygote count in true genotype: {hetero_count}")
    print(f"Homo to hetero ratio: {homo_count / hetero_count}")


def main():
    calc_switch_error_rate()


if __name__ == "__main__":
    main()
