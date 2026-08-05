"""Module for comparing accuracy reports and summarizing numeric changes."""

import csv
import math
import sys
import os
from collections import OrderedDict

from src.accuracy_core import ACCURACY_REPORT_VALUE_NAMES
from src.accuracy_visualization import load_accuracy_report


DEFAULT_BASELINE_REPORT = "accu_report.txt"
DEFAULT_CURRENT_REPORT = os.path.join(
    "tests", "accuracy_tests", "reports_benchmark", "accu_report.txt"
)
RUNTIME_METRIC = "elapsed_seconds"


def _record_key(record):
    """Return the identity used to match records across reports."""

    return (record["file_name"], record["label"], record["metric_name"])


def _index_records(records):
    """Index report records by file, label, and metric."""

    indexed = OrderedDict()
    for record in records:
        if record["metric_name"] == RUNTIME_METRIC:
            continue
        indexed[_record_key(record)] = record

    return indexed


def _value_name(index):
    """Return the report value name for a value-list index."""

    if index < len(ACCURACY_REPORT_VALUE_NAMES):
        return ACCURACY_REPORT_VALUE_NAMES[index]

    return f"value_{index}"


def _is_comparable_number(value):
    """Return whether a value can be used in numeric comparisons."""

    return isinstance(value, float) and math.isfinite(value)


def _new_metric_summary(metric_name):
    """Create an empty summary accumulator for one metric."""

    return {
        "metric_name": metric_name,
        "matched_records": 0,
        "compared_values": 0,
        "unchanged_values": 0,
        "increased_values": 0,
        "decreased_values": 0,
        "skipped_values": 0,
        "baseline_sum": 0.0,
        "current_sum": 0.0,
        "delta_sum": 0.0,
        "abs_delta_sum": 0.0,
        "max_abs_delta": 0.0,
        "max_abs_delta_location": "",
    }


def _update_metric_summary(summary, baseline_record, current_record, tolerance):
    """Accumulate comparison statistics for one matched record."""

    summary["matched_records"] += 1
    values_to_compare = min(
        len(baseline_record["values"]), len(current_record["values"])
    )
    summary["skipped_values"] += abs(
        len(baseline_record["values"]) - len(current_record["values"])
    )

    for index in range(values_to_compare):
        baseline_value = baseline_record["values"][index]
        current_value = current_record["values"][index]

        if not (
            _is_comparable_number(baseline_value)
            and _is_comparable_number(current_value)
        ):
            summary["skipped_values"] += 1
            continue

        delta = current_value - baseline_value
        abs_delta = abs(delta)
        summary["compared_values"] += 1
        summary["baseline_sum"] += baseline_value
        summary["current_sum"] += current_value
        summary["delta_sum"] += delta
        summary["abs_delta_sum"] += abs_delta

        if abs_delta <= tolerance:
            summary["unchanged_values"] += 1
        elif delta > 0:
            summary["increased_values"] += 1
        else:
            summary["decreased_values"] += 1

        if abs_delta > summary["max_abs_delta"]:
            summary["max_abs_delta"] = abs_delta
            summary["max_abs_delta_location"] = ":".join(
                [
                    baseline_record["file_name"],
                    baseline_record["label"],
                    _value_name(index),
                ]
            )


def compare_accuracy_reports(
    baseline_report_path=DEFAULT_BASELINE_REPORT,
    current_report_path=DEFAULT_CURRENT_REPORT,
    tolerance=1e-12,
):
    """Compare two accuracy reports and summarize numeric changes by metric."""

    baseline_records = _index_records(load_accuracy_report(baseline_report_path))
    current_records = _index_records(load_accuracy_report(current_report_path))

    summaries = OrderedDict()
    baseline_keys = set(baseline_records)
    current_keys = set(current_records)
    matched_keys = sorted(baseline_keys & current_keys)

    for key in matched_keys:
        metric_name = key[2]
        if metric_name not in summaries:
            summaries[metric_name] = _new_metric_summary(metric_name)

        _update_metric_summary(
            summaries[metric_name],
            baseline_records[key],
            current_records[key],
            tolerance,
        )

    return {
        "summaries": list(summaries.values()),
        "missing_from_current": sorted(baseline_keys - current_keys),
        "missing_from_baseline": sorted(current_keys - baseline_keys),
    }


def _format_float(value):
    """Format floats compactly for terminal and CSV output."""

    return f"{value:.10g}"


def _finalize_summary(summary):
    """Return a printable summary row with derived averages."""

    compared_values = summary["compared_values"]
    if compared_values == 0:
        mean_baseline = math.nan
        mean_current = math.nan
        mean_delta = math.nan
        mean_abs_delta = math.nan
    else:
        mean_baseline = summary["baseline_sum"] / compared_values
        mean_current = summary["current_sum"] / compared_values
        mean_delta = summary["delta_sum"] / compared_values
        mean_abs_delta = summary["abs_delta_sum"] / compared_values

    return {
        "metric_name": summary["metric_name"],
        "matched_records": summary["matched_records"],
        "compared_values": compared_values,
        "unchanged_values": summary["unchanged_values"],
        "increased_values": summary["increased_values"],
        "decreased_values": summary["decreased_values"],
        "skipped_values": summary["skipped_values"],
        "mean_baseline": mean_baseline,
        "mean_current": mean_current,
        "mean_delta": mean_delta,
        "mean_abs_delta": mean_abs_delta,
        "max_abs_delta": summary["max_abs_delta"],
        "max_abs_delta_location": summary["max_abs_delta_location"],
    }


def _write_csv(rows, output_path):
    """Write finalized summary rows to CSV."""

    fieldnames = [
        "metric_name",
        "matched_records",
        "compared_values",
        "unchanged_values",
        "increased_values",
        "decreased_values",
        "skipped_values",
        "mean_baseline",
        "mean_current",
        "mean_delta",
        "mean_abs_delta",
        "max_abs_delta",
        "max_abs_delta_location",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _print_summary(rows, missing_from_current, missing_from_baseline, output_file):
    """Print finalized summary rows."""

    print(
        "metric_name,matched_records,compared_values,mean_delta,"
        "mean_abs_delta,max_abs_delta,max_abs_delta_location",
        file=output_file,
    )
    for row in rows:
        print(
            ",".join(
                [
                    row["metric_name"],
                    str(row["matched_records"]),
                    str(row["compared_values"]),
                    _format_float(row["mean_delta"]),
                    _format_float(row["mean_abs_delta"]),
                    _format_float(row["max_abs_delta"]),
                    row["max_abs_delta_location"],
                ]
            ),
            file=output_file,
        )

    if missing_from_current:
        print(
            f"\nRecords missing from current report: {len(missing_from_current)}",
            file=output_file,
        )
    if missing_from_baseline:
        print(
            f"Records missing from baseline report: {len(missing_from_baseline)}",
            file=output_file,
        )


def compare_reports(
    baseline_report_path=DEFAULT_BASELINE_REPORT,
    current_report_path=DEFAULT_CURRENT_REPORT,
    output_path=None,
    tolerance=1e-12,
    output_file=None,
):
    """Compare report metrics using function arguments."""

    comparison = compare_accuracy_reports(
        baseline_report_path=baseline_report_path,
        current_report_path=current_report_path,
        tolerance=tolerance,
    )
    rows = [_finalize_summary(summary) for summary in comparison["summaries"]]
    if output_file is None:
        output_file = sys.stdout

    _print_summary(
        rows,
        comparison["missing_from_current"],
        comparison["missing_from_baseline"],
        output_file,
    )

    if output_path is not None:
        _write_csv(rows, output_path)
