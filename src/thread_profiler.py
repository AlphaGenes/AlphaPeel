"""Profile AlphaPeel runtime across n_thread_fam and n_thread_loci settings.

The profiler runs AlphaPeel in fresh Python subprocesses for each thread-count
combination, writes raw and summary CSV files, and saves 3D runtime surfaces.
"""

import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from src.accuracy_core import SIMULATION_FIXTURE_DIR, generate_accuracy_argv


DEFAULT_THREAD_COUNTS = (1, 2, 4, 8)
DEFAULT_OUTPUT_ROOT = os.path.join("tests", "accuracy_tests", "outputs_thread_grid")
INNER_MODE = "inner"


def _normalize_thread_counts(thread_counts):
    """Return thread counts as a tuple of positive integers."""

    counts = tuple(int(count) for count in thread_counts)
    if not counts or any(count < 1 for count in counts):
        raise ValueError("thread_counts must contain positive integers")
    return counts


def set_arg(argv, option, value):
    """Set or append a two-token AlphaPeel CLI option."""

    argv = list(argv)
    if option in argv:
        argv[argv.index(option) + 1] = str(value)
    else:
        argv.extend([option, str(value)])
    return argv


def build_tinypeel_argv(
    sim_path,
    method,
    output_dir,
    n_thread_fam,
    n_thread_loci,
    n_cycle=5,
    seq_file=True,
    est_geno_error_prob=True,
    est_seq_error_prob=True,
    est_start_alt_allele_prob=False,
    est_alt_allele_prob=False,
):
    """Build the AlphaPeel argv for one profiling run."""

    argv = generate_accuracy_argv(
        sim_path,
        method,
        est_start_alt_allele_prob,
        est_geno_error_prob,
        est_seq_error_prob,
        seq_file,
        False,
        est_alt_allele_prob,
        False,
        False,
        str(output_dir),
    )
    argv = set_arg(argv, "-n_cycle", n_cycle)
    argv = set_arg(argv, "-n_thread_fam", n_thread_fam)
    argv = set_arg(argv, "-n_thread_loci", n_thread_loci)
    return argv


def run_one(
    sim_path,
    method,
    output_dir,
    n_thread_fam,
    n_thread_loci,
    n_cycle=5,
    seq_file=True,
    est_geno_error_prob=True,
    est_seq_error_prob=True,
    est_start_alt_allele_prob=False,
    est_alt_allele_prob=False,
    python_executable=None,
    working_dir=None,
):
    """Run one AlphaPeel subprocess and return elapsed seconds."""

    argv = build_tinypeel_argv(
        sim_path,
        method,
        output_dir,
        n_thread_fam,
        n_thread_loci,
        n_cycle=n_cycle,
        seq_file=seq_file,
        est_geno_error_prob=est_geno_error_prob,
        est_seq_error_prob=est_seq_error_prob,
        est_start_alt_allele_prob=est_start_alt_allele_prob,
        est_alt_allele_prob=est_alt_allele_prob,
    )
    if python_executable is None:
        python_executable = sys.executable
    if working_dir is None:
        working_dir = Path.cwd()

    command = [python_executable, "-m", "src.tinypeel.tinypeel", *argv]
    start = time.perf_counter()
    result = subprocess.run(
        command,
        cwd=working_dir,
        text=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        check=False,
    )
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        raise RuntimeError(
            "AlphaPeel profiling run failed\n"
            f"n_thread_fam={n_thread_fam} n_thread_loci={n_thread_loci}\n"
            f"exit code={result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
    return elapsed


def read_existing_keys(raw_csv):
    """Return completed run keys from an existing raw CSV."""

    if not raw_csv.exists():
        return set()

    keys = set()
    with raw_csv.open(newline="") as file_in:
        for row in csv.DictReader(file_in):
            keys.add(
                (
                    row["mode"],
                    int(row["n_thread_fam"]),
                    int(row["n_thread_loci"]),
                    int(row["replicate"]),
                )
            )
    return keys


def append_raw_row(raw_csv, row):
    """Append one raw timing row."""

    write_header = not raw_csv.exists()
    with raw_csv.open("a", newline="") as file_out:
        writer = csv.DictWriter(
            file_out,
            fieldnames=(
                "mode",
                "n_thread_fam",
                "n_thread_loci",
                "replicate",
                "runtime_seconds",
            ),
        )
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def load_raw_rows(raw_csv):
    """Load raw timing rows."""

    rows = []
    with raw_csv.open(newline="") as file_in:
        for row in csv.DictReader(file_in):
            rows.append(
                {
                    "mode": row["mode"],
                    "n_thread_fam": int(row["n_thread_fam"]),
                    "n_thread_loci": int(row["n_thread_loci"]),
                    "replicate": int(row["replicate"]),
                    "runtime_seconds": float(row["runtime_seconds"]),
                }
            )
    return rows


def write_summary(summary_csv, rows):
    """Write mean and standard deviation by grid cell."""

    groups = {}
    for row in rows:
        key = (row["mode"], row["n_thread_fam"], row["n_thread_loci"])
        groups.setdefault(key, []).append(row["runtime_seconds"])

    with summary_csv.open("w", newline="") as file_out:
        writer = csv.DictWriter(
            file_out,
            fieldnames=(
                "mode",
                "n_thread_fam",
                "n_thread_loci",
                "replicates",
                "runtime_mean_seconds",
                "runtime_std_seconds",
                "runtime_min_seconds",
            ),
        )
        writer.writeheader()
        for key in sorted(groups):
            values = np.array(groups[key], dtype=np.float64)
            writer.writerow(
                {
                    "mode": key[0],
                    "n_thread_fam": key[1],
                    "n_thread_loci": key[2],
                    "replicates": len(values),
                    "runtime_mean_seconds": float(np.mean(values)),
                    "runtime_std_seconds": float(np.std(values)),
                    "runtime_min_seconds": float(np.min(values)),
                }
            )


def summary_grid(rows, mode, thread_counts):
    """Return X, Y, Z grids for one mode."""

    means = {}
    for row in rows:
        if row["mode"] != mode:
            continue
        key = (row["n_thread_fam"], row["n_thread_loci"])
        means.setdefault(key, []).append(row["runtime_seconds"])

    x_values = np.array(thread_counts, dtype=np.float64)
    y_values = np.array(thread_counts, dtype=np.float64)
    x_grid, y_grid = np.meshgrid(x_values, y_values)
    z_grid = np.full_like(x_grid, np.nan, dtype=np.float64)
    for y_index, n_thread_loci in enumerate(thread_counts):
        for x_index, n_thread_fam in enumerate(thread_counts):
            values = means.get((n_thread_fam, n_thread_loci))
            if values:
                z_grid[y_index, x_index] = np.mean(values)
    return x_grid, y_grid, z_grid


def plot_surface(output_root, rows, thread_counts, mode=INNER_MODE):
    """Create a 3D runtime surface plot for one mode."""

    mpl_config_dir = output_root / "matplotlib_config"
    xdg_cache_dir = output_root / "cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    xdg_cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = output_root / f"thread_grid_runtime_surface_{mode}.png"
    figure = plt.figure(figsize=(7, 6))
    axis = figure.add_subplot(1, 1, 1, projection="3d")
    draw_surface(axis, rows, mode, thread_counts)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def draw_surface(axis, rows, mode, thread_counts):
    """Draw one runtime surface."""

    x_grid, y_grid, z_grid = summary_grid(rows, mode, thread_counts)
    axis.plot_surface(
        x_grid,
        y_grid,
        z_grid,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.4,
        alpha=0.9,
    )
    axis.scatter(x_grid, y_grid, z_grid, color="black", s=18)
    axis.set_title(f"{mode} locus parallel mode")
    axis.set_xlabel("n_thread_fam")
    axis.set_ylabel("n_thread_loci")
    axis.set_zlabel("runtime seconds")
    axis.set_xticks(thread_counts)
    axis.set_yticks(thread_counts)


def run_thread_profiler(
    thread_counts=DEFAULT_THREAD_COUNTS,
    replicates=1,
    n_cycle=5,
    output_root=DEFAULT_OUTPUT_ROOT,
    sim_path=SIMULATION_FIXTURE_DIR,
    method="multi",
    seq_file=True,
    est_geno_error_prob=True,
    est_seq_error_prob=True,
    est_start_alt_allele_prob=False,
    est_alt_allele_prob=False,
    keep_outputs=False,
    skip_existing=False,
    python_executable=None,
    working_dir=None,
):
    """Run the thread-count profiling grid and return generated artifacts."""

    thread_counts = _normalize_thread_counts(thread_counts)
    if replicates < 1:
        raise ValueError("replicates must be at least 1")
    if method != "multi":
        raise ValueError("thread profiling currently supports method='multi' only")

    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    raw_csv = output_root / "thread_grid_raw.csv"
    summary_csv = output_root / "thread_grid_summary.csv"
    completed = read_existing_keys(raw_csv) if skip_existing else set()

    mode = INNER_MODE
    for n_thread_fam in thread_counts:
        for n_thread_loci in thread_counts:
            for replicate in range(1, replicates + 1):
                key = (mode, n_thread_fam, n_thread_loci, replicate)
                if key in completed:
                    print(f"Skipping existing {key}")
                    continue
                run_dir = (
                    output_root
                    / "runs"
                    / mode
                    / f"n_fam{n_thread_fam}_n_loci{n_thread_loci}_rep{replicate}"
                )
                if run_dir.exists():
                    shutil.rmtree(run_dir)
                run_dir.mkdir(parents=True)
                print(
                    "Running "
                    f"n_thread_fam={n_thread_fam} "
                    f"n_thread_loci={n_thread_loci} replicate={replicate}"
                )
                elapsed = run_one(
                    sim_path,
                    method,
                    run_dir,
                    n_thread_fam,
                    n_thread_loci,
                    n_cycle=n_cycle,
                    seq_file=seq_file,
                    est_geno_error_prob=est_geno_error_prob,
                    est_seq_error_prob=est_seq_error_prob,
                    est_start_alt_allele_prob=est_start_alt_allele_prob,
                    est_alt_allele_prob=est_alt_allele_prob,
                    python_executable=python_executable,
                    working_dir=working_dir,
                )
                append_raw_row(
                    raw_csv,
                    {
                        "mode": mode,
                        "n_thread_fam": n_thread_fam,
                        "n_thread_loci": n_thread_loci,
                        "replicate": replicate,
                        "runtime_seconds": elapsed,
                    },
                )
                print(f"  runtime_seconds={elapsed:.6f}")
                if not keep_outputs:
                    shutil.rmtree(run_dir)

    rows = load_raw_rows(raw_csv)
    write_summary(summary_csv, rows)
    plot_path = plot_surface(output_root, rows, thread_counts, INNER_MODE)
    print(f"Wrote raw timings: {raw_csv}")
    print(f"Wrote summary: {summary_csv}")
    print(f"Wrote plot: {plot_path}")
