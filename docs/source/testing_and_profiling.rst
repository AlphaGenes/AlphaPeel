.. _testing-and-profiling:

Testing and profiling
=====================

Before committing code changes, run the tests/checks from a package built
from your modified source. 

This page collects the common checks used during
AlphaPeel development: the regular ``pytest`` suite,
accuracy benchmark, and profiling and coverage checks.

``pytest`` suite
~~~~~~~~~~~~~~~~

The ``pytest`` suite is ran automatically on the GitHub Actions CI for every pull request, 
but it is recommended to run the tests locally before committing code changes.

Install test dependencies
-------------------------
Install ``pytest`` and ``pytest-benchmark`` before running the ``pytest`` suite.
For more background, see the `pytest documentation
<https://docs.pytest.org/en/stable/getting-started.html>`_ and the
`pytest-benchmark documentation
<https://pytest-benchmark.readthedocs.io/en/latest/installation.html>`_.

.. code-block:: bash

    pip install pytest
    pip install pytest-benchmark

Build and reinstall the modified package
----------------------------------------

Run tests against the package built from your current source tree.
This is required for the ``tests/accuracy_tests/run_accu_test.py`` as 
it runs the tests on the built package rather than the source code.

.. code-block:: bash

    python -m build
    python -m pip uninstall AlphaPeel -y
    python -m pip install dist/*.whl

For more detail on building a local distribution, see :ref:`dist-install`.

Run the ``pytest`` suite
------------------------

From the repository root, run:

.. code-block:: bash

    pytest

If the tests pass, the beginning of the output should look similar to this:

.. code-block::

    ============================================================================================ test session starts =============================================================================================
    platform darwin -- Python 3.11.11, pytest-9.0.1, pluggy-1.6.0
    benchmark: 5.2.3 (defaults: timer=time.perf_counter disable_gc=False min_rounds=5 min_time=0.000005 max_time=1.0 calibration_precision=10 warmup=False warmup_iterations=100000)
    rootdir: /Users/xtang3/AlphaPeel
    configfile: pyproject.toml
    plugins: benchmark-5.2.3, memray-1.9.0
    collected 16 items

    tests/accuracy_tests/run_accu_test.py ...                                                                                                                                                              [ 18%]
    tests/functional_tests/run_func_test.py .............                                                                                                                                                  [100%]
    ...

Run a targeted functional test
------------------------------

To run one functional test, provide the pytest node id. For example, to run
``test_files``:

.. code-block:: bash

    pytest tests/functional_tests/run_func_test.py::TestClass::test_files

Inspect failed tests and printed output
---------------------------------------

When a test fails and you need to inspect printed output, add
``-s``:

.. code-block:: bash

    pytest -s

Run accuracy benchmarks
~~~~~~~~~~~~~~~~~~~~~~~

In addition to running ``pytest``, you can run the full accuracy benchmark suite. 
This is useful when you want to compare the impact of a
new change on accuracy and runtime across the configured benchmark cases in detail.

From the repository root, run:

.. code-block:: python

    from src.accuracy_runner import run_full_accuracy_suite

    run_full_accuracy_suite(run_name="benchmark")

This writes outputs under ``tests/accuracy_tests/outputs_benchmark`` and writes
the accuracy report to ``tests/accuracy_tests/reports_benchmark/accu_report.txt``.
The report is a comma-separated text file with records in the form
``file_name,label,metric_name,value``. For most accuracy metrics, ``value``
contains the population metric followed by the metrics for generations 1 to 5.
Runtime is written separately as ``runtime,<label>,elapsed_seconds,<seconds>``.

Here is an example of the first 3 lines of the accuracy report:

.. code-block::

    runtime,single,elapsed_seconds,10.043625167018035
    dosage,single,marker_corr,['0.6792', '0.5175', '0.8658', '0.8722', '0.869', '0.7148']
    dosage,single,ind_corr,['0.861', '0.6016', '0.9377', '0.9416', '0.9401', '0.8841']

The accuracy report would have a few hundred lines and could be difficult to read in a text editor. 
You can visualise the report with the
``create_accuracy_report_visualizations`` function described below.

It is recommended to save a copy of the accuracy report before the change and after the change, 
and compare them with the ``compare_reports`` function described below.

Visualise accuracy reports
--------------------------

You can create plots from an accuracy report:

.. code-block:: python

    from src.accuracy_visualization import create_accuracy_report_visualizations

    create_accuracy_report_visualizations(
        "tests/accuracy_tests/reports_benchmark/accu_report.txt",
        "tests/accuracy_tests/reports_benchmark/plots",
    )

By default, this creates plots for ``marker_corr``, ``abs_diff``, and
``correct_rate``, together with a separate runtime plot. To focus on a different
set of metrics or methods, pass explicit selections:

.. code-block:: python

    create_accuracy_report_visualizations(
        "tests/accuracy_tests/reports_benchmark/accu_report.txt",
        "tests/accuracy_tests/reports_benchmark/plots",
        metric_names=("marker_corr", "switch_error_rate", "phase_error_rate"),
        labels=("single", "multi", "hybrid"),
    )

These plots are intended as a development aid before opening a pull request.

Compare accuracy reports
------------------------

After you have a baseline report and a report from your modified code, 
you can use the ``compare_reports`` function to compare them. 
The comparison highlights changes in runtime and accuracy, 
and can show whether differences are concentrated in
particular generations.

You can compare two accuracy reports with the
``compare_reports`` function by running the following 
code from the repository root:

.. code-block:: python

    from src.accuracy_report_comparison import compare_reports
    compare_reports()

It by default compares the two reports ``accu_report.txt``
and ``tests/accuracy_tests/reports_benchmark/accu_report.txt``, but you 
can pass explicit paths to the two reports to compare with the arguments
``baseline_report_path`` and ``current_report_path``.

An example output of the comparison is:

.. code-block:: python

    metric_name,matched_records,compared_values,mean_delta,mean_abs_delta,max_abs_delta,max_abs_delta_location
    abs_diff,98,572,5.48548139e-06,5.664273926e-06,0.002599332529,hap_0.5:multi_est_alt_allele_prob:generation_1
    ind_corr,98,572,-6.468531469e-06,6.468531469e-06,0.0029,hap_0.5:multi_est_alt_allele_prob:generation_1
    marker_corr,98,572,-1.590909091e-05,1.835664336e-05,0.0081,hap_0.5:multi_est_alt_allele_prob:generation_1
    correct_heterozygote_rate,15,90,1.333333333e-07,2.666666667e-07,7.5e-06,hap_0.5:multi_est_alt_allele_prob:generation_2
    heterozygote_count,15,90,0,0,0,
    homo_to_hetero_ratio,15,90,0,0,0,
    homozygote_count,15,90,0,0,0,
    phase_error_rate,15,90,1.73e-05,1.736666667e-05,0.0013,hap_0.5:multi_est_alt_allele_prob:generation_1
    switch_error_rate,15,90,-1.900950475e-06,1.900950475e-06,0.0001425712856,hap_0.5:multi_est_alt_allele_prob:generation_1
    uncalled_rate,15,90,2e-07,5.222222222e-07,1.75e-05,hap_0.5:multi_est_start_alt_allele_prob_est_geno_error_prob_est_seq_error_prob_seq_file:generation_1
    wrong_homozygote_rate,15,90,-3.666666667e-07,4.888888889e-07,1.75e-05,hap_0.5:multi_est_start_alt_allele_prob_est_geno_error_prob_est_seq_error_prob_seq_file:generation_1
    correct_rate,8,32,7.063383194e-09,7.89150818e-09,1.850000031e-08,seg_prob:multi_est_start_alt_allele_prob_est_geno_error_prob_est_seq_error_prob:generation_5

You can pass the ``output_path`` argument to save the comparison results to a CSV file.

Profiling and coverage
~~~~~~~~~~~~~~~~~~~~~~

Profiling is useful when you need to identify performance bottlenecks.

Memory profiling
----------------

For memory profiling, create a small temporary Python driver ``benchmark_memray.py``
that runs the workflow you want to inspect. For example:

.. code-block:: python

    from src.accuracy_runner import run_full_accuracy_suite

    run_full_accuracy_suite(run_name="benchmark")

Save the driver outside version control, or remove it after profiling. Then run
it with ``memray``:

.. code-block:: bash

    pip install memray
    memray run benchmark_memray.py

This creates a ``.bin`` memory profile file. Inspect it by generating a
flamegraph:

.. code-block:: bash

    memray flamegraph memray-*.bin

This generates an HTML file that you can open in a browser to inspect memory usage.

Runtime profiling
-----------------

For line-wise runtime profiling, create a small temporary Python driver with
``line_profiler``. Add the functions you want to inspect, then run the accuracy
case inside the profiler:

.. code-block:: python

    from line_profiler import LineProfiler

    from src.accuracy_core import get_params, sim_path
    from src.accuracy_runner import run_accuracy_case
    from src.tinypeel.tinypeel import main, peeling_cycle, run_peeling_cycles

    lp = LineProfiler()

    for fn in [
        run_peeling_cycles,
        peeling_cycle,
        main,
    ]:
        lp.add_function(fn)

    # This is a single accuracy case from the benchmark suite. You can change it to
    # any other case you want to profile.
    case = (
        "multi",  # method
        False,    # est_start_alt_allele_prob
        True,     # est_geno_error_prob
        True,     # est_seq_error_prob
        True,     # seq_file
        False,    # alt_allele_prob_file
        False,    # est_alt_allele_prob
        False,    # metafounder
        False,    # x_chr
    )

    with lp:
        run_accuracy_case(
            get_params(),
            sim_path(),
            *case,
            benchmark=None,
            run_name="profile_line",
            TINYPEEL_DIRECT_IS_WARMED_UP=True,
        )

    lp.print_stats()

Save this driver outside version control, or remove it after profiling. If the
temporary file is named ``line_profile.py``, run:

.. code-block:: bash

    pip install line_profiler
    python line_profile.py

The example above skips a warmup run by setting
``TINYPEEL_DIRECT_IS_WARMED_UP=True``. In a fresh Python process, this means the
reported line timings include JIT compilation time. If you have
``TINYPEEL_DIRECT_IS_WARMED_UP=False``, then the reported line timings include one warmup run, 
and the actual runtime of the function, which is equivalent to the compilation time 
plus two times the runtime.

To keep a reference profile for later comparison, redirect the output to a text
file with a descriptive name:

.. code-block:: bash

    python line_profile.py > before_change_line_profile.txt

.. _multi-threading-benchmarking:

Multi-threading benchmarking
----------------------------

Use ``src.thread_profiler`` to benchmark the ``multi`` method across
``n_thread_fam`` and ``n_thread_loci`` combinations. The profiler runs each
grid cell in a fresh Python subprocess, records raw timings, writes a summary
CSV, and creates 3D runtime surface plots.

From the repository root, run:

.. code-block:: python

    from src.thread_profiler import run_thread_profiler

    run_thread_profiler()

By default, this benchmarks thread counts ``1, 2, 4, 8`` with one replicate per
grid cell and writes artifacts under
``tests/accuracy_tests/outputs_thread_grid``. To use a smaller or larger grid,
pass explicit arguments:

.. code-block:: python

    run_thread_profiler(
        thread_counts=(1, 2, 4, 8, 16),
        replicates=3,
        n_cycle=5,
        skip_existing=True,
    )

The raw timings are written to ``thread_grid_raw.csv``, the aggregated timings
to ``thread_grid_summary.csv``, and the runtime surface plot to
``thread_grid_runtime_surface_inner.png`` in the output directory. Keep the CSV
files from a previous run if you want a reference for later comparisons.
Because each grid cell is run in a fresh Python subprocess, timings include
process startup and first-run JIT compilation time. Use the same settings when
comparing results across code changes.

Test coverage
-------------

Coverage can be collected with ``coverage``:

.. code-block:: bash

    pip install coverage
    coverage run -m pytest tests/functional_tests
    coverage report

To generate an HTML coverage report, run:

.. code-block:: bash

    coverage html

.. note::

    Coverage reports cannot identify calls to JIT-compiled functions.
