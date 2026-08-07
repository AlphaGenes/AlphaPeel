"""Module for running accuracy benchmarks and assessments of AlphaPeel."""

from dataclasses import dataclass
import os

import numpy as np

from src.accuracy_assessment import (
    AssessmentContext,
    _write_accuracy_metric,
    assess_accuracy,
    assess_test_accuracy,
)
from src.accuracy_core import (
    AccuracyCase,
    AccuracyInputOptions,
    AccuracyRunPaths,
    benchmark_tinypeel_direct,
    build_accuracy_case_name,
    build_accuracy_output_path,
    build_accuracy_report_path,
    _fixture_file_path,
    generate_accuracy_argv,
    generate_command,
    get_accuracy_benchmark_cases,
    get_accuracy_benchmark_output_root,
    get_accuracy_benchmark_report_root,
    get_params,
    prepare_directory,
    run_command,
    run_tinypeel_direct,
    sim_path as default_sim_path,
    ensure_directory,
    warmup_tinypeel_direct,
)


@dataclass(frozen=True)
class AccuracyRunnerOptions:
    """Execution options for one accuracy case run."""

    benchmark: object = None
    run_name: str = "test_accu"
    tinypeel_direct_is_warmed_up: bool = False


def _read_map_marker_names(path):
    """Return marker names from an AlphaPeel map file."""

    with open(path, "r", encoding="utf-8") as file:
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
    simulation_path,
    seq_file,
    metafounder,
    x_chr,
    output_path,
):
    """Create subset inputs for the first stage of a hybrid benchmark."""

    ensure_directory(_hybrid_subset_input_dir(output_path))

    map_path = _fixture_file_path(
        simulation_path,
        "map_file",
        metafounder=metafounder,
        x_chr=x_chr,
    )
    seg_map_path = _fixture_file_path(
        simulation_path,
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
        simulation_path,
        locus_file_name,
        metafounder=metafounder,
        x_chr=x_chr,
    )
    subset_path = os.path.join(
        _hybrid_subset_input_dir(output_path),
        f"{locus_file_name}.txt",
    )
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


def _case_for_method(case, method):
    """Return ``case`` options with a different AlphaPeel method."""

    return AccuracyCase(
        method,
        est_start_alt_allele_prob=case.est_start_alt_allele_prob,
        est_geno_error_prob=case.est_geno_error_prob,
        est_seq_error_prob=case.est_seq_error_prob,
        seq_file=case.seq_file,
        alt_allele_prob_file=case.alt_allele_prob_file,
        est_alt_allele_prob=case.est_alt_allele_prob,
        metafounder=case.metafounder,
        x_chr=case.x_chr,
    )


def run_accuracy_case(  # pylint: disable=too-many-locals
    parameters,
    simulation_path,
    case,
    options=None,
):
    """Run AlphaPeel and evaluate outputs based on the specified parameters."""

    if options is None:
        options = AccuracyRunnerOptions()

    name = build_accuracy_case_name(case)
    output_path = build_accuracy_output_path(name, options.run_name)
    paths = AccuracyRunPaths(simulation_path=simulation_path, output_path=output_path)
    prepare_directory(output_path)
    runtime_seconds = None

    if case.method == "hybrid":
        multi_output_path = _hybrid_multi_output_path(output_path)
        multi_paths = AccuracyRunPaths(
            simulation_path=simulation_path,
            output_path=multi_output_path,
        )
        prepare_directory(multi_output_path)
        multi_file_overrides = _prepare_hybrid_multi_inputs(
            simulation_path,
            case.seq_file,
            case.metafounder,
            case.x_chr,
            output_path,
        )
        hybrid_file_overrides = {
            "seg_file": os.path.join(multi_output_path, ".seg_prob.txt")
        }
        multi_options = AccuracyInputOptions(
            file_overrides=multi_file_overrides,
            extra_input_files=["map_file"],
        )
        hybrid_options = AccuracyInputOptions(
            file_overrides=hybrid_file_overrides,
            input_files_method="hybrid",
        )

        if options.run_name == "test_accu":
            commands = [
                generate_command(
                    multi_paths,
                    _case_for_method(case, "multi"),
                    input_options=multi_options,
                ),
                generate_command(
                    paths,
                    _case_for_method(case, "single"),
                    input_options=hybrid_options,
                ),
            ]
            options.benchmark(_run_commands, commands)
        else:
            argvs = [
                generate_accuracy_argv(
                    multi_paths,
                    _case_for_method(case, "multi"),
                    input_options=multi_options,
                ),
                generate_accuracy_argv(
                    paths,
                    _case_for_method(case, "single"),
                    input_options=hybrid_options,
                ),
            ]
            _run_tinypeel_direct_sequence(argvs)

    elif options.run_name == "test_accu":
        command = generate_command(
            paths,
            case,
        )
        options.benchmark(run_command, command)

    else:
        argv = generate_accuracy_argv(
            paths,
            case,
        )

        if not options.tinypeel_direct_is_warmed_up:
            warmup_tinypeel_direct(argv, output_path=output_path)

        runtime_seconds = benchmark_tinypeel_direct(argv)

    report_path = build_accuracy_report_path(options.run_name)
    ensure_directory(os.path.dirname(report_path))
    with open(report_path, "a", encoding="utf-8") as file_out:
        if options.run_name == "test_accu":
            assess_test_accuracy(
                simulation_path,
                parameters,
                output_path,
                case.method,
                file_out,
            )
        elif options.run_name == "benchmark":
            if runtime_seconds is not None:
                _write_accuracy_metric(
                    file_out,
                    "runtime",
                    name,
                    "elapsed_seconds",
                    runtime_seconds,
                )

            assess_accuracy(
                AssessmentContext(
                    simulation_path=simulation_path,
                    parameters=parameters,
                    output_path=output_path,
                    name=name,
                    method=case.method,
                    file_out=file_out,
                    metafounder=case.metafounder,
                    x_chr=case.x_chr,
                )
            )


def run_full_accuracy_suite(run_name="benchmark"):
    """Run the full direct-call accuracy benchmark suite."""

    tinypeel_direct_is_warmed_up = False

    prepare_directory(get_accuracy_benchmark_output_root(run_name))
    prepare_directory(get_accuracy_benchmark_report_root(run_name))

    for case in get_accuracy_benchmark_cases():
        run_accuracy_case(
            get_params(),
            default_sim_path(),
            case,
            options=AccuracyRunnerOptions(
                benchmark=None,
                run_name=run_name,
                tinypeel_direct_is_warmed_up=tinypeel_direct_is_warmed_up,
            ),
        )
        if not tinypeel_direct_is_warmed_up:
            tinypeel_direct_is_warmed_up = True
