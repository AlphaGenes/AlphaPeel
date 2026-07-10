import ast
import os

import numpy as np

from src.accuracy_core import (
    ACCURACY_REPORT_VALUE_NAMES,
    DEFAULT_VISUALIZATION_METRICS,
    ensure_directory,
)


def _parse_report_float(value):
    """Parse numeric report values, preserving nan as numpy.nan."""

    if isinstance(value, str) and value.lower() == "nan":
        return np.nan

    return float(value)


def _parse_accuracy_report_value(value):
    """Parse a scalar or list-valued accuracy report field."""

    value = value.strip()
    if value.startswith("["):
        return [_parse_report_float(item) for item in ast.literal_eval(value)]

    return [_parse_report_float(value)]


def load_accuracy_report(report_path):
    """Load an accuracy report into a list of structured records.

    The report uses four comma-separated fields, but list-valued metric fields
    also contain commas. Splitting each line at the first three commas preserves
    the current report format without requiring quoted CSV output.
    """

    records = []
    with open(report_path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line:
                continue

            fields = line.split(",", 3)
            if len(fields) != 4:
                raise ValueError(
                    f"Expected 4 fields at {report_path}:{line_number}, "
                    f"found {len(fields)}."
                )

            file_name, label, metric_name, raw_value = fields
            values = _parse_accuracy_report_value(raw_value)
            record = {
                "file_name": file_name,
                "label": label,
                "metric_name": metric_name,
                "values": values,
                "population": values[0],
            }
            for index, value_name in enumerate(ACCURACY_REPORT_VALUE_NAMES):
                if index < len(values):
                    record[value_name] = values[index]

            records.append(record)

    return records


def _matches_optional_filter(value, allowed_values):
    """Return whether a value passes an optional allow-list filter."""

    return allowed_values is None or value in allowed_values


def filter_accuracy_records(
    records,
    metric_names=DEFAULT_VISUALIZATION_METRICS,
    file_names=None,
    labels=None,
):
    """Filter loaded accuracy records for visualization."""

    metric_names = set(metric_names) if metric_names is not None else None
    file_names = set(file_names) if file_names is not None else None
    labels = set(labels) if labels is not None else None

    return [
        record
        for record in records
        if _matches_optional_filter(record["metric_name"], metric_names)
        and _matches_optional_filter(record["file_name"], file_names)
        and _matches_optional_filter(record["label"], labels)
    ]


def _unique_in_order(values):
    """Return unique values while preserving first-seen order."""

    unique_values = []
    seen = set()
    for value in values:
        if value not in seen:
            unique_values.append(value)
            seen.add(value)

    return unique_values


def _safe_plot_name(*parts):
    """Return a filesystem-friendly plot filename stem."""

    return "_".join(
        "".join(char if char.isalnum() else "_" for char in str(part)).strip("_")
        for part in parts
        if part
    )


def _get_pyplot():
    """Import matplotlib lazily for report visualization helpers."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "matplotlib is required for accuracy report visualizations."
        ) from error

    return plt


def _plot_output_path(output_dir, *name_parts):
    """Build an output path for a generated visualization."""

    ensure_directory(output_dir)
    return os.path.join(output_dir, f"{_safe_plot_name(*name_parts)}.png")


def _plot_colors(plt, count):
    """Return a distinct color sequence for a plot."""

    if count <= 0:
        return []

    color_map_name = "tab20" if count <= 20 else "nipy_spectral"
    color_map = plt.get_cmap(color_map_name, count)
    return [color_map(index) for index in range(count)]


def plot_accuracy_population_heatmaps(
    report_path,
    output_dir,
    metric_names=DEFAULT_VISUALIZATION_METRICS,
    file_names=None,
    labels=None,
):
    """Create one population-level heatmap per selected metric."""

    plt = _get_pyplot()
    records = filter_accuracy_records(
        load_accuracy_report(report_path),
        metric_names=metric_names,
        file_names=file_names,
        labels=labels,
    )
    generated_paths = []

    for metric_name in _unique_in_order(record["metric_name"] for record in records):
        metric_records = [
            record for record in records if record["metric_name"] == metric_name
        ]
        metric_labels = _unique_in_order(record["label"] for record in metric_records)
        metric_files = _unique_in_order(
            record["file_name"] for record in metric_records
        )
        values = {
            (record["label"], record["file_name"]): record["population"]
            for record in metric_records
        }
        matrix = np.array(
            [
                [values.get((label, file_name), np.nan) for file_name in metric_files]
                for label in metric_labels
            ],
            dtype=float,
        )

        max_label_length = max(len(label) for label in metric_labels)
        figure_width = max(10, len(metric_files) * 1.4 + max_label_length * 0.12)
        figure_height = max(4.5, len(metric_labels) * 0.35)
        figure, axis = plt.subplots(
            figsize=(figure_width, figure_height),
            constrained_layout=True,
        )
        image = axis.imshow(matrix, aspect="auto", cmap="viridis")
        axis.set_title(f"{metric_name} population metric")
        axis.set_xlabel("file_name")
        axis.set_ylabel("label")
        axis.set_xticks(np.arange(len(metric_files)))
        axis.set_xticklabels(metric_files, rotation=45, ha="right")
        axis.set_yticks(np.arange(len(metric_labels)))
        axis.set_yticklabels(metric_labels)
        figure.colorbar(image, ax=axis)

        output_path = _plot_output_path(output_dir, "population_heatmap", metric_name)
        figure.savefig(output_path, dpi=200)
        plt.close(figure)
        generated_paths.append(output_path)

    return generated_paths


def plot_accuracy_generation_profiles(
    report_path,
    output_dir,
    metric_names=DEFAULT_VISUALIZATION_METRICS,
    file_names=None,
    labels=None,
):
    """Create generation profile line plots for selected metrics and files."""

    plt = _get_pyplot()
    records = filter_accuracy_records(
        load_accuracy_report(report_path),
        metric_names=metric_names,
        file_names=file_names,
        labels=labels,
    )
    generated_paths = []

    metric_file_pairs = _unique_in_order(
        (record["metric_name"], record["file_name"])
        for record in records
        if len(record["values"]) > 1
    )
    for metric_name, file_name in metric_file_pairs:
        plot_records = [
            record
            for record in records
            if record["metric_name"] == metric_name
            and record["file_name"] == file_name
            and len(record["values"]) > 1
        ]
        if not plot_records:
            continue

        plot_labels = [record["label"] for record in plot_records]
        max_label_length = max(len(label) for label in plot_labels)
        figure_width = min(18, max(10, 6 + max_label_length * 0.12))
        figure, axis = plt.subplots(figsize=(figure_width, 5.5))
        generations = np.arange(1, len(ACCURACY_REPORT_VALUE_NAMES))
        colors = _plot_colors(plt, len(plot_records))
        for index, record in enumerate(plot_records):
            generation_values = record["values"][1:]
            axis.plot(
                generations,
                generation_values,
                color=colors[index],
                marker="o",
                linestyle="dashed",
                label=record["label"],
            )

        axis.set_title(f"{metric_name} by generation: {file_name}")
        axis.set_xlabel("generation")
        axis.set_ylabel(metric_name)
        axis.set_xticks(generations)
        axis.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize="small")
        figure.tight_layout()

        output_path = _plot_output_path(
            output_dir,
            "generation_profile",
            metric_name,
            file_name,
        )
        figure.savefig(output_path, dpi=200)
        plt.close(figure)
        generated_paths.append(output_path)

    return generated_paths


def plot_accuracy_runtime(report_path, output_dir, labels=None):
    """Create a runtime bar plot for benchmark cases."""

    plt = _get_pyplot()
    labels = set(labels) if labels is not None else None
    records = [
        record
        for record in load_accuracy_report(report_path)
        if record["file_name"] == "runtime"
        and record["metric_name"] == "elapsed_seconds"
        and _matches_optional_filter(record["label"], labels)
    ]
    if not records:
        return []

    records = sorted(records, key=lambda record: record["population"])
    runtime_labels = [record["label"] for record in records]
    runtime_values = [record["population"] for record in records]
    max_label_length = max(len(label) for label in runtime_labels)
    figure_width = min(18, max(10, 6 + max_label_length * 0.12))
    figure_height = max(4.5, len(records) * 0.35)

    figure, axis = plt.subplots(
        figsize=(figure_width, figure_height),
        constrained_layout=True,
    )
    colors = _plot_colors(plt, len(records))
    positions = np.arange(len(records))
    axis.barh(positions, runtime_values, color=colors)
    axis.set_title("Runtime by benchmark case")
    axis.set_xlabel("elapsed_seconds")
    axis.set_ylabel("label")
    axis.set_yticks(positions)
    axis.set_yticklabels(runtime_labels)
    axis.invert_yaxis()

    output_path = _plot_output_path(output_dir, "runtime", "elapsed_seconds")
    figure.savefig(output_path, dpi=200)
    plt.close(figure)

    return [output_path]


def create_accuracy_report_visualizations(
    report_path,
    output_dir,
    metric_names=DEFAULT_VISUALIZATION_METRICS,
    file_names=None,
    labels=None,
    include_generation_profiles=True,
    include_runtime=True,
):
    """Create the default accuracy report visualization set."""

    generated_paths = plot_accuracy_population_heatmaps(
        report_path,
        output_dir,
        metric_names=metric_names,
        file_names=file_names,
        labels=labels,
    )

    if include_generation_profiles:
        generated_paths.extend(
            plot_accuracy_generation_profiles(
                report_path,
                output_dir,
                metric_names=metric_names,
                file_names=file_names,
                labels=labels,
            )
        )

    if include_runtime:
        generated_paths.extend(
            plot_accuracy_runtime(
                report_path,
                output_dir,
                labels=labels,
            )
        )

    return generated_paths
