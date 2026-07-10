import pytest
import os
import shutil
import numpy as np
from src.utils import (
    get_accuracy_benchmark_output_root,
    get_accuracy_benchmark_report_root,
    get_params,
    prepare_directory,
)


@pytest.fixture(scope="session")
def test_get_params():
    return get_params()


def pytest_configure(config):
    """
    Prepare path and report file for accuracy tests
    """

    prepare_directory(get_accuracy_benchmark_output_root("test_accu"))
    prepare_directory(get_accuracy_benchmark_report_root("test_accu"))


@pytest.hookimpl()
def pytest_terminal_summary(terminalreporter):
    try:
        data = np.genfromtxt(
            os.path.join(
                "tests", "accuracy_tests", "reports_test_accu", "accu_report.txt"
            ),
            delimiter=",",
            names=["file", "method", "metric", "value"],
            dtype=[
                ("file", "U20"),
                ("method", "U20"),
                ("metric", "U30"),
                ("value", float),
            ],
        )
    except FileNotFoundError:
        return

    terminalreporter.write_sep("=", " Accuracy")

    files = np.unique(data["file_"])
    metrics = np.unique(data["metric"])

    for metric in metrics:
        terminalreporter.write_sep("~", metric)
        metric_data = data[data["metric"] == metric]

        if metric == "abs_diff":
            terminalreporter.write_line(
                "The metric abs_diff is the sum of absolute difference divided by the sum of the number of loci being counted. "
            )
            terminalreporter.write_line("The lower the value, the better the accuracy.")
        elif metric == "marker_corr":
            terminalreporter.write_line(
                "Pearson correlation evaluated at markers, ranged from 0 to 1."
            )
        elif metric == "ind_corr":
            terminalreporter.write_line(
                "Pearson correlation evaluated at individuals, ranged from 0 to 1."
            )
        elif metric == "correct_rate":
            terminalreporter.write_line(
                "Summing up the probabilities of the true state from the output data divided by the number of loci being counted."
            )
        elif metric == "switch_error_rate":
            terminalreporter.write_line(
                "Calculation follows the definition of SER from: https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X"
            )
        elif metric == "phase_error_rate":
            terminalreporter.write_line(
                "Calculation follows the definition of PER_intra from: https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X"
            )
        elif metric == "uncalled_rate":
            terminalreporter.write_line(
                "The number of loci with haplotypes uncalled divided by the total number of loci."
            )
        elif metric == "wrong_homozygote_rate":
            terminalreporter.write_line(
                "The number of loci that are called homozygote but are actually heterozygote, divided by the total number of loci."
            )
        elif metric == "true heterozygote rate":
            terminalreporter.write_line(
                "The number of loci that are indeed heterozygote and imputed heterozygote, divided by the total number of loci. This is approximately the maximum possible number of switch error rate."
            )

        bar_char = "#"
        empty_char = "."

        for file in files:
            file_data = metric_data[metric_data["file_"] == file]
            if len(file_data) != 0:
                terminalreporter.write_sep("-", file)
                terminalreporter.write_line("{:<20} {:<10}".format("Method", "Value"))
                terminalreporter.write_sep("-")

                values = file_data["value"]
                max_value = np.max(values)
                min_value = np.min(values)
                cols, _ = shutil.get_terminal_size()
                max_length = int(cols * 0.7)

                for row in file_data:
                    value = row["value"]
                    if max_value == min_value:
                        bar_length = max_length
                    else:
                        bar_length = int(
                            max_length * (value - min_value) / (max_value - min_value)
                        )
                    bar = bar_char * bar_length + empty_char * (max_length - bar_length)
                    terminalreporter.write_line(
                        "{:<20} {:.3f} | {} |".format(row["method"], value, bar)
                    )
