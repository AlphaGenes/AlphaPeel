import os

import numpy as np

from src.accuracy_assessment import (
    _write_accuracy_metric,
    assess_accuracy,
    assess_test_accuracy,
)
from src.accuracy_core import (
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
    sim_path,
    ensure_directory,
    warmup_tinypeel_direct,
)


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
    TINYPEEL_DIRECT_IS_WARMED_UP=False,
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
    runtime_seconds = None

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

        if not TINYPEEL_DIRECT_IS_WARMED_UP:
            warmup_tinypeel_direct(argv, output_path=output_path)

        runtime_seconds = benchmark_tinypeel_direct(argv)

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
            if runtime_seconds is not None:
                _write_accuracy_metric(
                    file_out,
                    "runtime",
                    name,
                    "elapsed_seconds",
                    runtime_seconds,
                )

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

    TINYPEEL_DIRECT_IS_WARMED_UP = False

    prepare_directory(get_accuracy_benchmark_output_root(run_name))
    prepare_directory(get_accuracy_benchmark_report_root(run_name))

    for case in get_accuracy_benchmark_cases():
        run_accuracy_case(
            get_params(),
            sim_path(),
            *case,
            benchmark=None,
            run_name=run_name,
            TINYPEEL_DIRECT_IS_WARMED_UP=TINYPEEL_DIRECT_IS_WARMED_UP,
        )
        if not TINYPEEL_DIRECT_IS_WARMED_UP:
            TINYPEEL_DIRECT_IS_WARMED_UP = True
