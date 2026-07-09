import argparse
import os
import shutil
import subprocess
import warnings

import numpy as np

import src.tinypeel.tinypeel as tinypeel


ACCURACY_TEST_ROOT = os.path.join("tests", "accuracy_tests")
SIMULATION_PARAMETERS_FILE = os.path.join(
    ACCURACY_TEST_ROOT, "simulation_parameters.txt"
)
SIMULATION_FIXTURE_DIR = os.path.join(ACCURACY_TEST_ROOT, "sim_for_alphapeel_accu_test")

DEFAULT_ALPHA_PEEL_ARGS = {
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

BASE_ACCURACY_FILES = [
    "dosage",
    "geno_0.333",
    "hap_0.5",
    "geno_prob",
    "phased_geno_prob",
]
METAFOUNDER_ACCURACY_FILES = [
    "dosage",
    "geno_prob",
    "phased_geno_prob",
]
ROWS_PER_INDIVIDUAL = {
    "dosage": 1,
    "geno_0.333": 1,
    "hap_0.5": 2,
    "geno_prob": 3,
    "phased_geno_prob": 4,
    "seg_prob": 4,
}
SEG_PROB_START_GEN = 2


def get_params():
    with open(SIMULATION_PARAMETERS_FILE, "r") as file:
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

    return SIMULATION_FIXTURE_DIR


def get_accuracy_benchmark_output_root(run_name="test_accu"):
    """Return the root directory for direct-call accuracy benchmark outputs."""

    return os.path.join(ACCURACY_TEST_ROOT, f"outputs_{run_name}")


def get_accuracy_benchmark_report_root(run_name="test_accu"):
    """Return the root directory for direct-call accuracy benchmark reports."""

    return os.path.join(ACCURACY_TEST_ROOT, f"reports_{run_name}")


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
    use_shell = isinstance(command, str)
    result = subprocess.run(
        command,
        shell=use_shell,
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
        str(param)
        for param in [
            method,
            est_start_alt_allele_prob,
            est_geno_error_prob,
            est_seq_error_prob,
            seq_file,
            alt_allele_prob_file,
            est_alt_allele_prob,
            metafounder,
            x_chr,
        ]
        if param
    )


def _alpha_peel_arguments(method):
    """Return the common AlphaPeel arguments for an accuracy run."""

    return {
        "method": method,
        **DEFAULT_ALPHA_PEEL_ARGS,
    }


def _input_files_for_case(method, seq_file=False, alt_allele_prob_file=False):
    """Return the fixture input file keys needed for an accuracy case."""

    input_files = ["ped_file", "seq_file" if seq_file else "geno_file"]

    if alt_allele_prob_file:
        input_files.append("alt_allele_prob_file")
    if method == "hybrid":
        input_files.extend(["map_file", "seg_map_file", "seg_file"])

    return input_files


def _fixture_file_prefix(metafounder=False, x_chr=False):
    """Return the fixture filename prefix for special benchmark modes."""

    if metafounder:
        return "metafounder_"
    if x_chr:
        return "X_chr_"
    return ""


def _fixture_file_path(sim_path, file_name, metafounder=False, x_chr=False):
    """Build the path to one benchmark fixture file."""

    prefix = _fixture_file_prefix(metafounder=metafounder, x_chr=x_chr)
    return os.path.join(sim_path, f"{prefix}{file_name}.txt")


def _input_file_path(
    sim_path,
    file_name,
    metafounder=False,
    x_chr=False,
    file_overrides=None,
):
    """Build an input file path, allowing a case to override selected files."""

    if file_overrides and file_name in file_overrides:
        return file_overrides[file_name]

    return _fixture_file_path(
        sim_path,
        file_name,
        metafounder=metafounder,
        x_chr=x_chr,
    )


def _add_argument(argv, key, value):
    """Append a CLI argument to ``argv`` when building a direct-call run."""

    argv.append(f"-{key}")
    if value is not None:
        argv.append(value)


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
    file_overrides=None,
    extra_input_files=None,
    input_files_method=None,
):
    """Generate an ``argv`` list for a direct ``tinypeel.main`` accuracy run."""

    argv = []
    input_files = _input_files_for_case(
        input_files_method or method,
        seq_file,
        alt_allele_prob_file,
    )
    if extra_input_files:
        input_files.extend(extra_input_files)

    arguments = _alpha_peel_arguments(method)

    if est_start_alt_allele_prob:
        arguments["est_start_alt_allele_prob"] = None
    if est_geno_error_prob and est_seq_error_prob:
        arguments["est_geno_error_prob"] = None
        arguments["est_seq_error_prob"] = None
    if est_alt_allele_prob:
        arguments["est_alt_allele_prob"] = None
    for file_name in input_files:
        argv.extend(
            [
                f"-{file_name}",
                _input_file_path(
                    sim_path,
                    file_name,
                    metafounder=metafounder,
                    x_chr=x_chr,
                    file_overrides=file_overrides,
                ),
            ]
        )

    if x_chr:
        argv.append("-x_chr")

    for key, value in arguments.items():
        _add_argument(argv, key, value)

    argv.extend(["-out_file", f"{output_path}{os.sep}"])

    return argv


def generate_command(
    sim_path,
    method,
    output_path,
    est_start_alt_allele_prob=False,
    est_geno_error_prob=False,
    est_seq_error_prob=False,
    seq_file=False,
    alt_allele_prob_file=False,
    est_alt_allele_prob=False,
    metafounder=False,
    x_chr=False,
    file_overrides=None,
    extra_input_files=None,
    input_files_method=None,
):
    """Generate the command used to run AlphaPeel for a test case.

    :param sim_path: Directory containing simulation input files.
    :type sim_path: str
    :param method: AlphaPeel method to run, such as ``single``, ``multi`` or
        ``hybrid``.
    :type method: str
    :param output_path: Output directory for AlphaPeel results.
    :type output_path: str
    :return: The command argv to pass to ``subprocess.run``.
    :rtype: list[str]
    """

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
        file_overrides=file_overrides,
        extra_input_files=extra_input_files,
        input_files_method=input_files_method,
    )

    return ["AlphaPeel", *argv]


def run_tinypeel_direct(argv):
    """Run AlphaPeel by calling ``src.tinypeel.tinypeel.main`` directly."""

    tinypeel.main(argv=argv)


def _accuracy_files(method, metafounder=False):
    """Return output file stems to compare for an accuracy run."""

    file_names = METAFOUNDER_ACCURACY_FILES if metafounder else BASE_ACCURACY_FILES
    file_names = list(file_names)

    if method == "multi":
        file_names.append("seg_prob")

    return file_names


def _accuracy_dimensions(get_params):
    """Return generation, individual, and locus dimensions from parameters."""

    n_gen = int(get_params["nGen"])
    n_ind_per_gen = int(get_params["nInd"] / n_gen)
    n_loci_all = int(get_params["nLociAll"])

    return n_gen, n_ind_per_gen, n_loci_all


def _load_accuracy_matrix(path, n_loci_all):
    """Load an AlphaPeel output or truth matrix for accuracy comparison."""

    return np.loadtxt(path, usecols=np.arange(1, n_loci_all + 1))


def _truth_file_path(sim_path, file_name, metafounder=False, x_chr=False):
    """Build the path to a truth file for an accuracy comparison."""

    if metafounder:
        true_file = f"true-metafounder_{file_name}.txt"
    elif x_chr:
        true_file = f"true-X_chr_{file_name}.txt"
    else:
        true_file = f"true-{file_name}.txt"

    return os.path.join(sim_path, true_file)


def _comparison_slice(matrix, file_name, n_ind_per_gen, n_row_per_ind):
    """Return the portion of a matrix used by the test accuracy report."""

    if file_name == "seg_prob":
        start = SEG_PROB_START_GEN * (n_ind_per_gen * n_row_per_ind)
        return matrix[start:, 1:]

    return matrix[:, 1:]


def _generation_slice(matrix, gen, n_ind_per_gen, n_row_per_ind):
    """Return rows belonging to one generation."""

    start = gen * (n_ind_per_gen * n_row_per_ind)
    end = (gen + 1) * (n_ind_per_gen * n_row_per_ind)

    return matrix[start:end]


def assess_test_accuracy(
    sim_path,
    get_params,
    output_path,
    method,
    file_out,
):
    file_to_check = _accuracy_files(method)
    _, n_ind_per_gen, n_loci_all = _accuracy_dimensions(get_params)

    for file_name in file_to_check:
        n_row_per_ind = ROWS_PER_INDIVIDUAL[file_name]
        file_path = os.path.join(output_path, f".{file_name}.txt")
        true_path = os.path.join(sim_path, f"true-{file_name}.txt")

        try:
            new_file = _load_accuracy_matrix(file_path, n_loci_all)
        except ValueError:
            print(f"Error loading {file_path}")
            continue

        try:
            true_file = _load_accuracy_matrix(true_path, n_loci_all)
        except ValueError:
            print(f"Error loading {true_path}")
            continue

        marker_corr = get_marker_corr(
            _comparison_slice(new_file, file_name, n_ind_per_gen, n_row_per_ind),
            _comparison_slice(true_file, file_name, n_ind_per_gen, n_row_per_ind),
        )

        file_out.write(f"{file_name},{method},marker_corr,{marker_corr}\n")

        if file_name == "seg_prob":
            ind_corr = get_ind_corr(
                new_file[:, :],
                true_file[:, :],
                n_ind_per_gen,
                n_row_per_ind,
                start_gen=SEG_PROB_START_GEN,
            )
        else:
            ind_corr = get_ind_corr(
                new_file[:, :],
                true_file[:, :],
                n_ind_per_gen,
                n_row_per_ind,
            )

        file_out.write(f"{file_name},{method},ind_corr,{ind_corr}\n")

        abs_diff = get_abs_diff(
            _comparison_slice(new_file, file_name, n_ind_per_gen, n_row_per_ind),
            _comparison_slice(true_file, file_name, n_ind_per_gen, n_row_per_ind),
            n_row_per_ind,
        )

        file_out.write(f"{file_name},{method},abs_diff,{abs_diff}\n")


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

    file_to_check = _accuracy_files(method, metafounder)
    n_gen, n_ind_per_gen, n_loci_all = _accuracy_dimensions(get_params)

    print(" ")
    print(f"Test: {name}")

    for file_name in file_to_check:
        n_row_per_ind = ROWS_PER_INDIVIDUAL[file_name]
        file_path = os.path.join(output_path, f".{file_name}.txt")
        true_path = _truth_file_path(
            sim_path,
            file_name,
            metafounder=metafounder,
            x_chr=x_chr,
        )

        new_file = _load_accuracy_matrix(file_path, n_loci_all)
        true_file = _load_accuracy_matrix(true_path, n_loci_all)

        if x_chr and file_name == "hap_0.5":
            new_file[true_file == 9] = 0
            true_file[true_file == 9] = 0

        marker_corr = [
            str(
                get_marker_corr(
                    _comparison_slice(
                        new_file, file_name, n_ind_per_gen, n_row_per_ind
                    ),
                    _comparison_slice(
                        true_file, file_name, n_ind_per_gen, n_row_per_ind
                    ),
                )
            )
        ]
        for gen in range(n_gen):
            if gen in [0, 1] and file_name == "seg_prob":
                marker_corr.append("nan")
                continue
            marker_corr.append(
                str(
                    get_marker_corr(
                        _generation_slice(
                            new_file,
                            gen,
                            n_ind_per_gen,
                            n_row_per_ind,
                        ),
                        _generation_slice(
                            true_file,
                            gen,
                            n_ind_per_gen,
                            n_row_per_ind,
                        ),
                    )
                )
            )

        file_out.write(f"{file_name},{name},marker_corr,{marker_corr}\n")

        if file_name == "seg_prob":
            ind_corr = [
                str(
                    get_ind_corr(
                        new_file[:, :],
                        true_file[:, :],
                        n_ind_per_gen,
                        n_row_per_ind,
                        start_gen=SEG_PROB_START_GEN,
                    )
                )
            ]
        else:
            ind_corr = [
                str(
                    get_ind_corr(
                        new_file[:, :],
                        true_file[:, :],
                        n_ind_per_gen,
                        n_row_per_ind,
                    )
                )
            ]
        for gen in range(n_gen):
            if gen in [0, 1] and file_name == "seg_prob":
                ind_corr.append("nan")
                continue
            ind_corr.append(
                str(
                    get_ind_corr(
                        new_file[:, :],
                        true_file[:, :],
                        n_ind_per_gen,
                        n_row_per_ind,
                        gen,
                        gen + 1,
                    )
                )
            )

        file_out.write(f"{file_name},{name},ind_corr,{ind_corr}\n")


def _read_map_marker_names(path):
    """Return marker names from an AlphaPeel map file."""

    with open(path, "r") as file:
        return [line.split()[1] for line in file if line.strip()]


def _subset_locus_file(source_path, target_path, columns):
    """Write a genotype or sequence file containing only selected loci."""

    data = np.loadtxt(source_path, dtype=np.int64)
    np.savetxt(target_path, data[:, columns], fmt="%d")


def _hybrid_multi_output_path(output_path):
    """Return the first-stage multi-locus output directory for a hybrid case."""

    return os.path.join(output_path, "multi_stage")


def _hybrid_subset_input_dir(output_path):
    """Return the directory for first-stage subset genotype or sequence inputs."""

    return os.path.join(output_path, "subset_inputs")


def _prepare_hybrid_multi_inputs(
    sim_path,
    seq_file,
    metafounder,
    x_chr,
    output_path,
):
    """Create subset inputs for the first stage of a hybrid benchmark."""

    subset_input_dir = _hybrid_subset_input_dir(output_path)
    ensure_directory(subset_input_dir)

    map_path = _fixture_file_path(
        sim_path,
        "map_file",
        metafounder=metafounder,
        x_chr=x_chr,
    )
    seg_map_path = _fixture_file_path(
        sim_path,
        "seg_map_file",
        metafounder=metafounder,
        x_chr=x_chr,
    )
    full_markers = _read_map_marker_names(map_path)
    subset_markers = _read_map_marker_names(seg_map_path)
    marker_columns = {marker: index + 1 for index, marker in enumerate(full_markers)}
    missing_markers = [
        marker for marker in subset_markers if marker not in marker_columns
    ]
    if missing_markers:
        raise ValueError(
            "Hybrid seg_map_file contains markers that are not in map_file: "
            f"{missing_markers}"
        )

    columns = [0] + [marker_columns[marker] for marker in subset_markers]

    locus_file_name = "seq_file" if seq_file else "geno_file"
    source_path = _fixture_file_path(
        sim_path,
        locus_file_name,
        metafounder=metafounder,
        x_chr=x_chr,
    )
    subset_path = os.path.join(subset_input_dir, f"{locus_file_name}.txt")
    _subset_locus_file(source_path, subset_path, columns)

    return {
        locus_file_name: subset_path,
        "map_file": seg_map_path,
    }


def _run_commands(commands):
    """Run a sequence of shell commands."""

    for command in commands:
        run_command(command)


def _run_tinypeel_direct_sequence(argvs):
    """Run a sequence of direct tinypeel calls."""

    for argv in argvs:
        run_tinypeel_direct(argv)


def run_accuracy_case(
    get_params,
    sim_path,
    method,
    est_start_alt_allele_prob=False,
    est_geno_error_prob=False,
    est_seq_error_prob=False,
    seq_file=False,
    alt_allele_prob_file=False,
    est_alt_allele_prob=False,
    metafounder=False,
    x_chr=False,
    benchmark=None,
    run_name="test_accu",
):
    """Run AlphaPeel and evaluate outputs based on the specified parameters."""

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
        multi_output_path = _hybrid_multi_output_path(output_path)
        prepare_directory(multi_output_path)
        multi_file_overrides = _prepare_hybrid_multi_inputs(
            sim_path,
            seq_file,
            metafounder,
            x_chr,
            output_path,
        )
        hybrid_file_overrides = {
            "seg_file": os.path.join(multi_output_path, ".seg_prob.txt")
        }

        if run_name == "test_accu":
            commands = [
                generate_command(
                    sim_path,
                    "multi",
                    multi_output_path,
                    est_start_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    est_alt_allele_prob,
                    metafounder,
                    x_chr,
                    file_overrides=multi_file_overrides,
                    extra_input_files=["map_file"],
                ),
                generate_command(
                    sim_path,
                    "single",
                    output_path,
                    est_start_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    est_alt_allele_prob,
                    metafounder,
                    x_chr,
                    file_overrides=hybrid_file_overrides,
                    input_files_method="hybrid",
                ),
            ]
            benchmark(_run_commands, commands)
        else:
            argvs = [
                generate_accuracy_argv(
                    sim_path,
                    "multi",
                    est_start_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    est_alt_allele_prob,
                    metafounder,
                    x_chr,
                    multi_output_path,
                    file_overrides=multi_file_overrides,
                    extra_input_files=["map_file"],
                ),
                generate_accuracy_argv(
                    sim_path,
                    "single",
                    est_start_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    est_alt_allele_prob,
                    metafounder,
                    x_chr,
                    output_path,
                    file_overrides=hybrid_file_overrides,
                    input_files_method="hybrid",
                ),
            ]
            _run_tinypeel_direct_sequence(argvs)

    elif run_name == "test_accu":
        command = generate_command(
            sim_path,
            method,
            output_path,
            est_start_alt_allele_prob,
            est_geno_error_prob,
            est_seq_error_prob,
            seq_file,
            alt_allele_prob_file,
            est_alt_allele_prob,
            metafounder,
            x_chr,
        )
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
        if isinstance(start_gen, int):
            if not isinstance(end_gen, int):
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
