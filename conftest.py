"""Pytest configuration file for tests."""

import os
import shutil
import pytest
import numpy as np
from src.accuracy_core import (
    get_accuracy_benchmark_output_root,
    get_accuracy_benchmark_report_root,
    get_params,
    prepare_directory,
)

ACCURACY_REPORT_PATH = os.path.join(
    "tests",
    "accuracy_tests",
    "reports_test_accu",
    "accu_report.txt",
)

METRIC_DESCRIPTIONS = {
    "abs_diff": [
        "The metric abs_diff is the sum of absolute difference divided by "
        "the sum of the number of loci being counted. ",
        "The lower the value, the better the accuracy.",
    ],
    "marker_corr": [
        "Pearson correlation evaluated at markers, ranged from 0 to 1.",
    ],
    "ind_corr": [
        "Pearson correlation evaluated at individuals, ranged from 0 to 1.",
    ],
    "correct_rate": [
        "Summing up the probabilities of the true state from the output data "
        "divided by the number of loci being counted.",
    ],
    "switch_error_rate": [
        "Calculation follows the definition of SER from: "
        "https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X",
    ],
    "phase_error_rate": [
        "Calculation follows the definition of PER_intra from: "
        "https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X",
    ],
    "uncalled_rate": [
        "The number of loci with haplotypes uncalled divided by the total "
        "number of loci.",
    ],
    "wrong_homozygote_rate": [
        "The number of loci that are called homozygote but are actually "
        "heterozygote, divided by the total number of loci.",
    ],
    "true heterozygote rate": [
        "The number of loci that are indeed heterozygote and imputed "
        "heterozygote, divided by the total number of loci. This is "
        "approximately the maximum possible number of switch error rate.",
    ],
}


@pytest.fixture(scope="session")
def test_get_params():
    """Simulation data parameters for accuracy tests."""
    return get_params()


# pylint: disable=unused-argument
def pytest_configure(config):
    """
    Prepare path and report file for accuracy tests
    """

    prepare_directory(get_accuracy_benchmark_output_root("test_accu"))
    prepare_directory(get_accuracy_benchmark_report_root("test_accu"))


def _load_accuracy_report():
    """Load the generated accuracy benchmark report."""

    return np.genfromtxt(
        ACCURACY_REPORT_PATH,
        delimiter=",",
        names=["file", "method", "metric", "value"],
        dtype=[
            ("file", "U20"),
            ("method", "U20"),
            ("metric", "U30"),
            ("value", float),
        ],
    )


def _write_metric_description(terminalreporter, metric):
    """Write explanatory text for one accuracy metric."""

    for description in METRIC_DESCRIPTIONS.get(metric, []):
        terminalreporter.write_line(description)


def _bar_length(value, min_value, max_value, max_length):
    """Return the scaled bar length for one value."""

    if max_value == min_value:
        return max_length
    return int(max_length * (value - min_value) / (max_value - min_value))


def _write_accuracy_rows(terminalreporter, file_data):
    """Write all accuracy rows for one file."""

    bar_char = "#"
    empty_char = "."
    values = file_data["value"]
    max_value = np.max(values)
    min_value = np.min(values)
    cols, _ = shutil.get_terminal_size()
    max_length = int(cols * 0.7)

    terminalreporter.write_line(f"{'Method':<20} {'Value':<10}")
    terminalreporter.write_sep("-")

    for row in file_data:
        value = row["value"]
        bar_length = _bar_length(value, min_value, max_value, max_length)
        out_bar = bar_char * bar_length + empty_char * (max_length - bar_length)
        terminalreporter.write_line(f"{row['method']:<20} {value:.3f} | {out_bar} |")


def _write_metric_summary(terminalreporter, metric, metric_data, files):
    """Write the accuracy summary for one metric."""

    terminalreporter.write_sep("~", metric)
    _write_metric_description(terminalreporter, metric)

    for file in files:
        file_data = metric_data[metric_data["file_"] == file]
        if len(file_data) == 0:
            continue
        terminalreporter.write_sep("-", file)
        _write_accuracy_rows(terminalreporter, file_data)


@pytest.hookimpl()
def pytest_terminal_summary(terminalreporter):
    """
    Generate a summary of the accuracy test results.

    :param terminalreporter: The terminal reporter object used to write the summary to the terminal.
    :type terminalreporter: pytest.terminal.TerminalReporter
    """
    try:
        data = _load_accuracy_report()
    except FileNotFoundError:
        return

    terminalreporter.write_sep("=", " Accuracy")

    files = np.unique(data["file_"])
    metrics = np.unique(data["metric"])

    for metric in metrics:
        metric_data = data[data["metric"] == metric]
        _write_metric_summary(terminalreporter, metric, metric_data, files)
