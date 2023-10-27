import pytest
import operator
import os
import shutil
import numpy as np


@pytest.fixture(scope="session")
def get_params():
    param_file = os.path.join("tests", "accuracy_tests", "simulation_parameters.txt")
    with open(param_file, "r") as file:
        sim_params = [line.strip().split() for line in file]

    params = {}
    for param_name, param_value in sim_params:
        params[param_name] = float(param_value)

    return params


def pytest_configure(config):
    """
    Prepare path and report file for accuracy tests
    """
    accu_output_path = os.path.join("tests", "accuracy_tests", "outputs")
    if os.path.exists(accu_output_path):
        shutil.rmtree(accu_output_path)
    os.mkdir(accu_output_path)

    report_path = os.path.join("tests", "accuracy_tests", "accu_report.txt")
    f = open(report_path, "w")
    f.close()


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport():
    out = yield
    report = out.get_result()
    if (
        report.nodeid[:37] == "tests/accuracy_tests/run_accu_test.py"
        and report.when == "call"
    ):
        stdout = report.sections
        if "multi" in report.nodeid:
            accu = stdout[-1][-1].split("\n")[-17:]
            num_file = 5
        else:
            accu = stdout[-1][-1].split("\n")[-14:]
            num_file = 4
        name = accu[0].split()[-1]
        with open("tests/accuracy_tests/accu_report.txt", "a") as file:
            for i in range(num_file):
                assessed_file = accu[i * 3 + 1].split()[-1]
                file.write(name + "|" + assessed_file + "$" + accu[i * 3 + 2] + "\n")
                file.write(name + "|" + assessed_file + "$" + accu[i * 3 + 3] + "\n")


@pytest.hookimpl()
def pytest_terminal_summary(terminalreporter):
    param_file = os.path.join("tests", "accuracy_tests", "simulation_parameters.txt")
    with open(param_file, "r") as file:
        sim_params = [line.strip().split() for line in file]

    params = {}
    for param_name, param_value in sim_params:
        params[param_name] = float(param_value)

    nGen = int(params["nGen"])

    file_types = ["dosages", "called.0.1", "called_phase.0.1", "haps", "seg"]
    columns = (
        "Type",
        "Population Accu",
        "Gen1 Accu",
        "Gen2 Accu",
        "Gen3 Accu",
        "Gen4 Accu",
        "Gen5 Accu",
    )
    dt = {"names": columns, "formats": ("S69", "f4", "f4", "f4", "f4", "f4", "f4")}
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

    mkr_accu_file = {}
    ind_accu_file = {}
    for file_type in file_types:
        mkr_accu_file[file_type] = list(
            filter(
                lambda x: (
                    x["Type"].decode("UTF-8").split("|")[1].split("$")[0] == file_type
                ),
                mkr_accu,
            )
        )
        ind_accu_file[file_type] = list(
            filter(
                lambda x: (
                    x["Type"].decode("UTF-8").split("|")[1].split("$")[0] == file_type
                ),
                ind_accu,
            )
        )

    sorted_mkr_accu_file = {}
    sorted_ind_accu_file = {}
    for file_type in file_types:
        sorted_mkr_accu_file[file_type] = sorted(
            mkr_accu_file[file_type],
            key=operator.itemgetter(*columns[1:]),
            reverse=True,
        )
        sorted_ind_accu_file[file_type] = sorted(
            ind_accu_file[file_type],
            key=operator.itemgetter(*columns[1:]),
            reverse=True,
        )

    format_first_row = "{:<35} " + "{:<20} " * (nGen + 1)
    format_row = "{:<35} " + "{:<20.3f} " * (nGen + 1)

    def write_report(accu_type, sorted=False):
        if not sorted:
            terminalreporter.write_sep("=", accu_type + " Accuaracy")
            if accu_type == "Marker":
                reports = mkr_accu_file
            else:
                reports = ind_accu_file
        else:
            terminalreporter.write_sep(
                "=", accu_type + " Accuaracy (Order by accuracy)"
            )
            if accu_type == "Marker":
                reports = sorted_mkr_accu_file
            else:
                reports = sorted_ind_accu_file

        for file_type in file_types:
            terminalreporter.write_sep("-", file_type)
            terminalreporter.write_line(format_first_row.format(*columns))
            for test in reports[file_type]:
                terminalreporter.write_line(
                    format_row.format(
                        test["Type"].decode("UTF-8").split("|")[0],
                        test["Population Accu"],
                        test["Gen1 Accu"],
                        test["Gen2 Accu"],
                        test["Gen3 Accu"],
                        test["Gen4 Accu"],
                        test["Gen5 Accu"],
                    )
                )

    write_report("Marker", sorted=False)
    write_report("Individual", sorted=False)
    write_report("Marker", sorted=True)
    write_report("Individual", sorted=True)
