"""Module for core accuracy assessment functionality."""

from dataclasses import dataclass
import os
import shutil
import subprocess
import time
from typing import Dict, List, Optional

from src.tinypeel import tinypeel


ACCURACY_TEST_ROOT = os.path.join("tests", "accuracy_tests")
SIMULATION_PARAMETERS_FILE = os.path.join(
    ACCURACY_TEST_ROOT, "simulation_parameters.txt"
)
SIMULATION_FIXTURE_DIR = os.path.join(ACCURACY_TEST_ROOT, "sim_for_alphapeel_accu_test")

DEFAULT_ALPHA_PEEL_ARGS = {
    "n_cycle": "5",
    "n_thread_fam": "2",
    "n_thread_loci": "3",
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
HAP_FILE = "hap_0.5"
DEFAULT_VISUALIZATION_METRICS = ("marker_corr", "abs_diff", "correct_rate")
ACCURACY_REPORT_VALUE_NAMES = (
    "population",
    "generation_1",
    "generation_2",
    "generation_3",
    "generation_4",
    "generation_5",
)


@dataclass(frozen=True)
class AccuracyCase:  # pylint: disable=too-many-instance-attributes
    """Configuration for one accuracy benchmark case."""

    method: str
    est_start_alt_allele_prob: bool = False
    est_geno_error_prob: bool = False
    est_seq_error_prob: bool = False
    seq_file: bool = False
    alt_allele_prob_file: bool = False
    est_alt_allele_prob: bool = False
    metafounder: bool = False
    x_chr: bool = False


@dataclass(frozen=True)
class AccuracyRunPaths:
    """Input and output paths for one accuracy run."""

    simulation_path: str
    output_path: str


@dataclass(frozen=True)
class AccuracyInputOptions:
    """Input-file overrides used when building AlphaPeel argv."""

    file_overrides: Optional[Dict[str, str]] = None
    extra_input_files: Optional[List[str]] = None
    input_files_method: Optional[str] = None


def get_params():
    """Get the simulation parameters from the simulation_parameters.txt file."""
    with open(SIMULATION_PARAMETERS_FILE, "r", encoding="utf-8") as file:
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
    """Run a command using subprocess.run and raise an error if it fails."""
    use_shell = isinstance(command, str)
    result = subprocess.run(
        command,
        shell=use_shell,
        text=True,
        capture_output=True,
        check=False,
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
        AccuracyCase("single"),
        AccuracyCase("single", est_start_alt_allele_prob=True),
        AccuracyCase("single", est_alt_allele_prob=True),
        AccuracyCase("multi"),
        AccuracyCase("multi", est_start_alt_allele_prob=True),
        AccuracyCase("multi", est_alt_allele_prob=True),
        AccuracyCase(
            "multi",
            est_start_alt_allele_prob=True,
            est_geno_error_prob=True,
            est_seq_error_prob=True,
        ),
        AccuracyCase("multi", seq_file=True),
        AccuracyCase("multi", est_start_alt_allele_prob=True, seq_file=True),
        AccuracyCase(
            "multi",
            est_start_alt_allele_prob=True,
            est_geno_error_prob=True,
            est_seq_error_prob=True,
            seq_file=True,
        ),
        AccuracyCase("hybrid"),
        AccuracyCase("hybrid", seq_file=True),
        AccuracyCase("single", alt_allele_prob_file=True, metafounder=True),
        AccuracyCase("single", est_start_alt_allele_prob=True, metafounder=True),
        AccuracyCase("single", est_alt_allele_prob=True, metafounder=True),
        AccuracyCase(
            "single",
            est_start_alt_allele_prob=True,
            est_alt_allele_prob=True,
            metafounder=True,
        ),
        AccuracyCase(
            "single",
            alt_allele_prob_file=True,
            est_alt_allele_prob=True,
            metafounder=True,
        ),
        AccuracyCase("single", x_chr=True),
        AccuracyCase("multi", x_chr=True),
        AccuracyCase("hybrid", x_chr=True),
    ]


def build_accuracy_case_name(case):
    """Build the benchmark case name from the accuracy case parameters."""

    return "_".join(
        str(param)
        for param in [
            case.method,
            "est_start_alt_allele_prob" if case.est_start_alt_allele_prob else None,
            "est_geno_error_prob" if case.est_geno_error_prob else None,
            "est_seq_error_prob" if case.est_seq_error_prob else None,
            "seq_file" if case.seq_file else None,
            "alt_allele_prob_file" if case.alt_allele_prob_file else None,
            "est_alt_allele_prob" if case.est_alt_allele_prob else None,
            "metafounder" if case.metafounder else None,
            "x_chr" if case.x_chr else None,
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


def _fixture_file_path(simulation_path, file_name, metafounder=False, x_chr=False):
    """Build the path to one benchmark fixture file."""

    prefix = _fixture_file_prefix(metafounder=metafounder, x_chr=x_chr)
    return os.path.join(simulation_path, f"{prefix}{file_name}.txt")


def _input_file_path(
    simulation_path,
    file_name,
    metafounder=False,
    x_chr=False,
    file_overrides=None,
):
    """Build an input file path, allowing a case to override selected files."""

    if file_overrides and file_name in file_overrides:
        return file_overrides[file_name]

    return _fixture_file_path(
        simulation_path,
        file_name,
        metafounder=metafounder,
        x_chr=x_chr,
    )


def _add_argument(argv, key, value):
    """Append a CLI argument to ``argv`` when building a direct-call run."""

    argv.append(f"-{key}")
    if value is not None:
        argv.append(value)


def generate_accuracy_argv(paths, case, input_options=None):
    """Generate an ``argv`` list for a direct ``tinypeel.main`` accuracy run."""

    if input_options is None:
        input_options = AccuracyInputOptions()

    argv = []
    input_files = _input_files_for_case(
        input_options.input_files_method or case.method,
        case.seq_file,
        case.alt_allele_prob_file,
    )
    if input_options.extra_input_files:
        input_files.extend(input_options.extra_input_files)

    arguments = _alpha_peel_arguments(case.method)

    if case.est_start_alt_allele_prob:
        arguments["est_start_alt_allele_prob"] = None
    if case.est_geno_error_prob and case.est_seq_error_prob:
        arguments["est_geno_error_prob"] = None
        arguments["est_seq_error_prob"] = None
    if case.est_alt_allele_prob:
        arguments["est_alt_allele_prob"] = None
    for file_name in input_files:
        argv.extend(
            [
                f"-{file_name}",
                _input_file_path(
                    paths.simulation_path,
                    file_name,
                    metafounder=case.metafounder,
                    x_chr=case.x_chr,
                    file_overrides=input_options.file_overrides,
                ),
            ]
        )

    if case.x_chr:
        argv.append("-x_chr")

    for key, value in arguments.items():
        _add_argument(argv, key, value)

    argv.extend(["-out_file", f"{paths.output_path}{os.sep}"])

    return argv


def generate_command(paths, case, input_options=None):
    """Generate the command used to run AlphaPeel for a test case.

    :param paths: Simulation input and output paths for the run.
    :type paths: AccuracyRunPaths
    :param case: AlphaPeel method and accuracy-case options.
    :type case: AccuracyCase
    :param input_options: Optional input-file overrides for special cases.
    :type input_options: AccuracyInputOptions or None
    :return: The command argv to pass to ``subprocess.run``.
    :rtype: list[str]
    """

    argv = generate_accuracy_argv(paths, case, input_options=input_options)

    return ["AlphaPeel", *argv]


def run_tinypeel_direct(argv):
    """Run AlphaPeel by calling ``src.tinypeel.tinypeel.main`` directly."""

    tinypeel.main(argv=argv)


def warmup_tinypeel_direct(argv, output_path=None):
    """Warm up JIT compilation."""

    run_tinypeel_direct(argv)
    if output_path is not None:
        prepare_directory(output_path)


def benchmark_tinypeel_direct(argv):
    """Run AlphaPeel and return elapsed seconds."""

    start_time = time.perf_counter()
    run_tinypeel_direct(argv)

    return time.perf_counter() - start_time
