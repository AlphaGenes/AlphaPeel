"""Module for assessing the accuracy of AlphaPeel output against truth files."""
from dataclasses import dataclass
import os
import warnings

import numpy as np

from src.accuracy_core import (
    BASE_ACCURACY_FILES,
    HAP_FILE,
    METAFOUNDER_ACCURACY_FILES,
    ROWS_PER_INDIVIDUAL,
    SEG_PROB_START_GEN,
)


@dataclass(frozen=True)
class AccuracyDimensions:
    """Simulation dimensions used by accuracy reports."""

    n_gen: int
    n_ind_per_gen: int
    n_loci_all: int


@dataclass(frozen=True)
class FileMetricContext:
    """Per-file dimensions used to calculate accuracy metrics."""

    file_name: str
    n_gen: int
    n_ind_per_gen: int
    n_row_per_ind: int


@dataclass(frozen=True)
class AssessmentContext:  # pylint: disable=too-many-instance-attributes
    """Inputs shared by one benchmark accuracy assessment."""

    simulation_path: str
    parameters: dict
    output_path: str
    name: str
    method: str
    file_out: object
    metafounder: bool
    x_chr: bool


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
    return AccuracyDimensions(
        n_gen=n_gen,
        n_ind_per_gen=int(get_params["nInd"] / n_gen),
        n_loci_all=int(get_params["nLociAll"]),
    )


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


def _load_accuracy_pair(output_path, sim_path, file_name, n_loci_all, **truth_kwargs):
    """Load one output/truth matrix pair for an accuracy comparison."""

    file_path = os.path.join(output_path, f".{file_name}.txt")
    true_path = _truth_file_path(sim_path, file_name, **truth_kwargs)

    return (
        _load_accuracy_matrix(file_path, n_loci_all),
        _load_accuracy_matrix(true_path, n_loci_all),
    )


def _test_truth_file_path(sim_path, file_name):
    """Build the truth path used by the pytest accuracy report."""

    return os.path.join(sim_path, f"true-{file_name}.txt")


def _load_test_accuracy_pair(output_path, sim_path, file_name, n_loci_all):
    """Load one output/truth matrix pair for the pytest accuracy report."""

    file_path = os.path.join(output_path, f".{file_name}.txt")
    true_path = _test_truth_file_path(sim_path, file_name)

    try:
        output = _load_accuracy_matrix(file_path, n_loci_all)
    except ValueError:
        print(f"Error loading {file_path}")
        return None

    try:
        truth = _load_accuracy_matrix(true_path, n_loci_all)
    except ValueError:
        print(f"Error loading {true_path}")
        return None

    return output, truth


def _mask_x_chr_hap_missing(output, truth, file_name, x_chr):
    """Mask missing X chromosome haplotype truth values before comparison."""

    if x_chr and file_name == HAP_FILE:
        output[truth == 9] = 0
        truth[truth == 9] = 0


def _normal_metric_output(output, file_name):
    """Return output values used by non-switch accuracy metrics."""

    if file_name != HAP_FILE:
        return output

    metric_output = output.astype(float, copy=True)
    metric_output[metric_output == 9] = np.nan

    return metric_output


def _overall_comparison(output, truth, file_name, n_ind_per_gen, n_row_per_ind):
    """Return output/truth slices used for whole-file metrics."""

    return (
        _comparison_slice(output, file_name, n_ind_per_gen, n_row_per_ind),
        _comparison_slice(truth, file_name, n_ind_per_gen, n_row_per_ind),
    )


def _ind_corr_for_file(output, truth, file_name, n_ind_per_gen, n_row_per_ind):
    """Return individual correlation with seg_prob generation handling."""

    if file_name == "seg_prob":
        return get_ind_corr(
            output,
            truth,
            n_ind_per_gen,
            n_row_per_ind,
            start_gen=SEG_PROB_START_GEN,
        )

    return get_ind_corr(output, truth, n_ind_per_gen, n_row_per_ind)


def _write_accuracy_metric(file_out, file_name, label, metric_name, value):
    """Write one metric line to an accuracy report."""

    file_out.write(f"{file_name},{label},{metric_name},{value}\n")


def _seg_prob_generation_is_skipped(file_name, gen):
    """Return whether a per-generation seg_prob metric should be reported as nan."""

    return file_name == "seg_prob" and gen < SEG_PROB_START_GEN


def _per_generation_metric(
    metric_func,
    output,
    truth,
    context,
    *metric_args,
):
    """Return a list of per-generation metric values as strings."""

    values = []
    for gen in range(context.n_gen):
        if _seg_prob_generation_is_skipped(context.file_name, gen):
            values.append("nan")
            continue

        values.append(
            str(
                metric_func(
                    _generation_slice(
                        output,
                        gen,
                        context.n_ind_per_gen,
                        context.n_row_per_ind,
                    ),
                    _generation_slice(
                        truth,
                        gen,
                        context.n_ind_per_gen,
                        context.n_row_per_ind,
                    ),
                    *metric_args,
                )
            )
        )

    return values


def _benchmark_metric_values(
    metric_func,
    output,
    truth,
    context,
    *metric_args,
):
    """Return whole-file and per-generation benchmark metric values."""

    overall_output, overall_truth = _overall_comparison(
        output,
        truth,
        context.file_name,
        context.n_ind_per_gen,
        context.n_row_per_ind,
    )
    overall_value = metric_func(overall_output, overall_truth, *metric_args)

    return [str(overall_value)] + _per_generation_metric(
        metric_func,
        output,
        truth,
        context,
        *metric_args,
    )


def _benchmark_ind_corr_values(output, truth, context):
    """Return whole-file and per-generation individual-correlation values."""

    overall_start_gen = None
    if context.file_name == "seg_prob":
        overall_start_gen = SEG_PROB_START_GEN

    values = [
        str(
            get_ind_corr(
                output,
                truth,
                context.n_ind_per_gen,
                context.n_row_per_ind,
                start_gen=overall_start_gen,
            )
        )
    ]

    return values + _per_generation_metric(
        get_ind_corr,
        output,
        truth,
        context,
        context.n_ind_per_gen,
        context.n_row_per_ind,
    )


def _safe_rate(numerator, denominator):
    """Return a rate, or nan when the denominator is empty."""

    if denominator == 0:
        return np.nan

    return numerator / denominator


def get_hap_switch_error_metrics(called_file, true_file, n_ind, n_loci_all):
    """Calculate switch, phase, and related haplotype call metrics."""

    counts = {
        "switch_error": 0,
        "phase_error": 0,
        "uncalled": 0,
        "wrong_homo": 0,
        "true_hetero": 0,
        "homo": 0,
        "hetero": 0,
    }

    for ind in range(n_ind):
        hap_p_new = called_file[ind * 2]
        hap_m_new = called_file[ind * 2 + 1]
        hap_p_true = true_file[ind * 2]
        hap_m_true = true_file[ind * 2 + 1]

        switched = False
        for loci in range(n_loci_all):
            if hap_p_true[loci] == hap_m_true[loci]:
                counts["homo"] += 1
            else:
                counts["hetero"] += 1

            if np.isnan(hap_p_new[loci]) or np.isnan(hap_m_new[loci]):
                counts["uncalled"] += 1
                continue

            if hap_p_new[loci] == 9 or hap_m_new[loci] == 9:
                counts["uncalled"] += 1
                continue

            if hap_p_true[loci] + hap_m_true[loci] == 1:
                if hap_p_new[loci] == hap_m_new[loci]:
                    counts["wrong_homo"] += 1
                    continue

                counts["true_hetero"] += 1
                if hap_p_new[loci] != hap_p_true[loci]:
                    counts["phase_error"] += 1

                if (hap_p_new[loci] != hap_p_true[loci] and switched is False) or (
                    hap_p_new[loci] == hap_p_true[loci] and switched is True
                ):
                    switched = not switched
                    counts["switch_error"] += 1

    genotype_count = n_ind * n_loci_all
    switch_opportunity_count = n_ind * (n_loci_all - 1)

    return [
        (
            "switch_error_rate",
            _safe_rate(counts["switch_error"], switch_opportunity_count),
        ),
        ("phase_error_rate", _safe_rate(counts["phase_error"], genotype_count)),
        ("uncalled_rate", _safe_rate(counts["uncalled"], genotype_count)),
        ("wrong_homozygote_rate", _safe_rate(counts["wrong_homo"], genotype_count)),
        (
            "correct_heterozygote_rate",
            _safe_rate(counts["true_hetero"], genotype_count),
        ),
        ("homozygote_count", counts["homo"]),
        ("heterozygote_count", counts["hetero"]),
        ("homo_to_hetero_ratio", _safe_rate(counts["homo"], counts["hetero"])),
    ]


def _benchmark_hap_switch_error_metrics(output, truth, n_gen, n_ind_per_gen):
    """Return whole-file and per-generation switch-error metrics."""

    n_loci_all = output.shape[1]
    metric_values = [
        (metric_name, [str(value)])
        for metric_name, value in get_hap_switch_error_metrics(
            output,
            truth,
            n_gen * n_ind_per_gen,
            n_loci_all,
        )
    ]
    values_by_metric = dict(metric_values)

    for gen in range(n_gen):
        gen_output = _generation_slice(
            output,
            gen,
            n_ind_per_gen,
            ROWS_PER_INDIVIDUAL[HAP_FILE],
        )
        gen_truth = _generation_slice(
            truth,
            gen,
            n_ind_per_gen,
            ROWS_PER_INDIVIDUAL[HAP_FILE],
        )

        for metric_name, value in get_hap_switch_error_metrics(
            gen_output,
            gen_truth,
            n_ind_per_gen,
            n_loci_all,
        ):
            values_by_metric[metric_name].append(str(value))

    return [
        (metric_name, values_by_metric[metric_name]) for metric_name, _ in metric_values
    ]


def _write_test_hap_switch_error_metrics(target, output, truth, dimensions):
    """Write test accuracy switch-error metrics for haplotype output."""

    file_out, file_name, method = target
    for metric_name, value in get_hap_switch_error_metrics(
        output,
        truth,
        dimensions.n_gen * dimensions.n_ind_per_gen,
        output.shape[1],
    ):
        if metric_name in [
            "homozygote_count",
            "heterozygote_count",
            "homo_to_hetero_ratio",
        ]:
            continue
        _write_accuracy_metric(file_out, file_name, method, metric_name, value)


def _write_benchmark_hap_switch_error_metrics(target, output, truth, dimensions):
    """Write benchmark switch-error metrics for haplotype output."""

    file_out, file_name, name = target
    for metric_name, values in _benchmark_hap_switch_error_metrics(
        output,
        truth,
        dimensions.n_gen,
        dimensions.n_ind_per_gen,
    ):
        _write_accuracy_metric(file_out, file_name, name, metric_name, values)


def _write_test_file_accuracy(target, dimensions, accuracy_pair):
    """Write test accuracy metrics for one output file."""

    file_out, file_name, method = target
    new_file, true_file = accuracy_pair
    n_row_per_ind = ROWS_PER_INDIVIDUAL[file_name]
    metric_new_file = _normal_metric_output(new_file, file_name)
    comparison_output, comparison_truth = _overall_comparison(
        metric_new_file,
        true_file,
        file_name,
        dimensions.n_ind_per_gen,
        n_row_per_ind,
    )
    _write_accuracy_metric(
        file_out,
        file_name,
        method,
        "marker_corr",
        get_marker_corr(comparison_output, comparison_truth),
    )
    _write_accuracy_metric(
        file_out,
        file_name,
        method,
        "ind_corr",
        _ind_corr_for_file(
            metric_new_file,
            true_file,
            file_name,
            dimensions.n_ind_per_gen,
            n_row_per_ind,
        ),
    )
    _write_accuracy_metric(
        file_out,
        file_name,
        method,
        "abs_diff",
        get_abs_diff(comparison_output, comparison_truth, n_row_per_ind),
    )

    if file_name == HAP_FILE:
        _write_test_hap_switch_error_metrics(
            target,
            new_file,
            true_file,
            dimensions,
        )
    if file_name == "seg_prob":
        correct_rate = get_correct_rate(comparison_output, comparison_truth)
        _write_accuracy_metric(
            file_out, file_name, method, "correct_rate", correct_rate
        )


def assess_test_accuracy(
    sim_path,
    get_params,
    output_path,
    method,
    file_out,
):
    """Assess pytest accuracy output against truth files."""

    dimensions = _accuracy_dimensions(get_params)
    for file_name in _accuracy_files(method):
        accuracy_pair = _load_test_accuracy_pair(
            output_path, sim_path, file_name, dimensions.n_loci_all
        )
        if accuracy_pair is None:
            continue
        _write_test_file_accuracy(
            (file_out, file_name, method),
            dimensions,
            accuracy_pair,
        )


def _benchmark_file_context(file_name, dimensions):
    """Build a file metric context from simulation dimensions."""

    return FileMetricContext(
        file_name=file_name,
        n_gen=dimensions.n_gen,
        n_ind_per_gen=dimensions.n_ind_per_gen,
        n_row_per_ind=ROWS_PER_INDIVIDUAL[file_name],
    )


def _write_benchmark_file_accuracy(context, dimensions, file_name):
    """Write benchmark accuracy metrics for one output file."""

    new_file, true_file = _load_accuracy_pair(
        context.output_path,
        context.simulation_path,
        file_name,
        dimensions.n_loci_all,
        metafounder=context.metafounder,
        x_chr=context.x_chr,
    )
    _mask_x_chr_hap_missing(new_file, true_file, file_name, context.x_chr)
    metric_new_file = _normal_metric_output(new_file, file_name)
    file_context = _benchmark_file_context(file_name, dimensions)

    _write_accuracy_metric(
        context.file_out,
        file_name,
        context.name,
        "marker_corr",
        _benchmark_metric_values(
            get_marker_corr,
            metric_new_file,
            true_file,
            file_context,
        ),
    )
    _write_accuracy_metric(
        context.file_out,
        file_name,
        context.name,
        "ind_corr",
        _benchmark_ind_corr_values(metric_new_file, true_file, file_context),
    )
    _write_accuracy_metric(
        context.file_out,
        file_name,
        context.name,
        "abs_diff",
        _benchmark_metric_values(
            get_abs_diff,
            metric_new_file,
            true_file,
            file_context,
            file_context.n_row_per_ind,
        ),
    )

    if file_name == HAP_FILE:
        _write_benchmark_hap_switch_error_metrics(
            (context.file_out, file_name, context.name),
            new_file,
            true_file,
            dimensions,
        )
    if file_name == "seg_prob":
        correct_rate = _benchmark_metric_values(
            get_correct_rate,
            new_file,
            true_file,
            file_context,
        )
        _write_accuracy_metric(
            context.file_out,
            file_name,
            context.name,
            "correct_rate",
            correct_rate,
        )


def assess_accuracy(context):
    """Assess benchmark output accuracy against truth files."""

    dimensions = _accuracy_dimensions(context.parameters)

    print(" ")
    print(f"Test: {context.name}")

    for file_name in _accuracy_files(context.method, context.metafounder):
        _write_benchmark_file_accuracy(context, dimensions, file_name)


# pylint: disable=too-many-arguments,too-many-positional-arguments
def assess_accuracy_from_parts(
    sim_path,
    get_params,
    output_path,
    name,
    method,
    file_out,
    metafounder,
    x_chr,
):
    """Backward-compatible wrapper around :func:`assess_accuracy`."""

    assess_accuracy(
        AssessmentContext(
            simulation_path=sim_path,
            parameters=get_params,
            output_path=output_path,
            name=name,
            method=method,
            file_out=file_out,
            metafounder=metafounder,
            x_chr=x_chr,
        )
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
        accus = []
        for i in range(real.shape[1]):
            valid = ~np.isnan(output[:, i]) & ~np.isnan(real[:, i])
            if np.sum(valid) < 2:
                accus.append(np.nan)
                continue

            accus.append(np.corrcoef(real[valid, i], output[valid, i])[0, 1])

        accus = np.array(accus)
        return round(np.nanmean(accus), 4)


# pylint: disable=invalid-name
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
        accus = []
        for i in range(real.shape[0]):
            valid = ~np.isnan(output[i, :]) & ~np.isnan(real[i, :])
            if np.sum(valid) < 2:
                accus.append(np.nan)
                continue

            accus.append(np.corrcoef(real[i, valid], output[i, valid])[0, 1])

        accus = np.array(accus)
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
    valid = ~np.isnan(output) & ~np.isnan(real)
    valid_count = np.sum(valid)
    if valid_count == 0:
        return np.nan

    return n_row_per_ind * np.sum(np.abs(output[valid] - real[valid])) / valid_count


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
