import pytest


@pytest.hookimpl(hookwrapper=True, tryfirst=True)
def pytest_runtest_makereport(item, call):
    out = yield
    report = out.get_result()

    if (
        report.when == "call"
        and report.nodeid[:37] == "tests/accuracy_tests/run_accu_test.py"
    ):
        print("---" * 35)
        print("nodeid: " + report.nodeid)
        accuracy_report = report.sections[-1][-1].split("\n")[-2]
        print(accuracy_report)
