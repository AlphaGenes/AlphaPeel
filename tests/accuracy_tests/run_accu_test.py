"""Accuracy tests for the AlphaPeel module."""

import pytest
from src.accuracy_core import AccuracyCase, sim_path
from src.accuracy_runner import AccuracyRunnerOptions, run_accuracy_case


@pytest.mark.parametrize(
    "method",
    [
        ("single"),
        ("multi"),
        ("hybrid"),
    ],
)
def test_accu(
    test_get_params,
    method,
    benchmark,
):
    """Run the accuracy benchmark across the configured input combinations.

    :param method: AlphaPeel method to run.
    :type method: str
    :param benchmark: Benchmark fixture used to time the command execution.
    :type benchmark: callable
    :return: None
    :rtype: None
    """

    run_accuracy_case(
        test_get_params,
        sim_path(),
        AccuracyCase(method),
        options=AccuracyRunnerOptions(
            benchmark=benchmark,
            run_name="test_accu",
        ),
    )
