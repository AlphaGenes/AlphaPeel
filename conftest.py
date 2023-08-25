import pytest
import operator
import numpy as np


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport():
    out = yield
    report = out.get_result()
    if (
        report.nodeid[:37] == "tests/accuracy_tests/run_accu_test.py"
        and report.when == "call"
    ):
        stdout = report.sections
        accu = stdout[-1][-1].split("\n")[-3:-1]
        with open("tests/accuracy_tests/accu_report.txt", "a") as file:
            file.write(report.nodeid[39:] + "$" + accu[0] + "\n")
            file.write(report.nodeid[39:] + "$" + accu[1] + "\n")


@pytest.hookimpl()
def pytest_terminal_summary(terminalreporter):
    columns = (
        "Type",
        "Population Accu",
        "Gen1 Accu",
        "Gen2 Accu",
        "Gen3 Accu",
        "Gen4 Accu",
        "Gen5 Accu",
    )
    dt = {"names": columns, "formats": ("S53", "f4", "f4", "f4", "f4", "f4", "f4")}
    accu = np.loadtxt("tests/accuracy_tests/accu_report.txt", dtype=dt)

    mkr_accu = list(
        filter(
            lambda x: (x["Type"].decode("UTF-8").split("$")[1] == "Marker_accuracies"),
            accu,
        )
    )
    ind_accu = list(
        filter(
            lambda x: (
                x["Type"].decode("UTF-8").split("$")[1] == "Individual_accuracies"
            ),
            accu,
        )
    )
    sorted_mkr_accu = sorted(
        mkr_accu, key=operator.itemgetter(*columns[1:]), reverse=True
    )
    sorted_ind_accu = sorted(
        ind_accu, key=operator.itemgetter(*columns[1:]), reverse=True
    )

    format_first_row = "{:<35} " + "{:<20} " * 6
    format_row = "{:<35} " + "{:<20.3f} " * 6

    terminalreporter.write_sep("-", "Marker Accuracy")
    terminalreporter.write_line(format_first_row.format(*columns))
    for test in mkr_accu:
        terminalreporter.write_line(
            format_row.format(
                test["Type"].decode("UTF-8").split("$")[0],
                test["Population Accu"],
                test["Gen1 Accu"],
                test["Gen2 Accu"],
                test["Gen3 Accu"],
                test["Gen4 Accu"],
                test["Gen5 Accu"],
            )
        )

    terminalreporter.write_sep("-", "Individual Accuracy")
    terminalreporter.write_line(format_first_row.format(*columns))
    for test in ind_accu:
        terminalreporter.write_line(
            format_row.format(
                test["Type"].decode("UTF-8").split("$")[0],
                test["Population Accu"],
                test["Gen1 Accu"],
                test["Gen2 Accu"],
                test["Gen3 Accu"],
                test["Gen4 Accu"],
                test["Gen5 Accu"],
            )
        )

    terminalreporter.write_sep("-", "Marker Accuracy (Order by accuracy)")
    terminalreporter.write_line(format_first_row.format(*columns))
    for test in sorted_mkr_accu:
        terminalreporter.write_line(
            format_row.format(
                test["Type"].decode("UTF-8").split("$")[0],
                test["Population Accu"],
                test["Gen1 Accu"],
                test["Gen2 Accu"],
                test["Gen3 Accu"],
                test["Gen4 Accu"],
                test["Gen5 Accu"],
            )
        )

    terminalreporter.write_sep("-", "Individual Accuracy (Order by accuracy)")
    terminalreporter.write_line(format_first_row.format(*columns))
    for test in sorted_ind_accu:
        terminalreporter.write_line(
            format_row.format(
                test["Type"].decode("UTF-8").split("$")[0],
                test["Population Accu"],
                test["Gen1 Accu"],
                test["Gen2 Accu"],
                test["Gen3 Accu"],
                test["Gen4 Accu"],
                test["Gen5 Accu"],
            )
        )
