#!/usr/bin/env python3
"""
This script was mostly written by Copilot with auto model.
"""

import io
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path


def load_matrix(path: str, delimiter: str = ',') -> np.ndarray:
    """Load a numeric matrix from a text file with a specified delimiter.
    The first column may be a string identifier and is skipped; remaining columns must be numeric.
    """
    with open(path) as fh:
        lines = fh.readlines()
    lines = lines[:-1]  # drop footer
    n_fields = len(lines[0].split(delimiter))
    data = np.loadtxt(io.StringIO("".join(lines)), delimiter=delimiter, usecols=range(1, n_fields))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def get_left_panel_points(root_dir: Path) -> list[tuple[float, float, str, str]]:
    """Collect (percent_h, percent_k, color, marker) points for all expected cohorts under root_dir."""
    # (color, marker, filled)
    style_map = {
        (10, 8): ("blue", "o", True),
        (10, 32): ("red", "o", True),
        (50, 8): ("blue", "o", False),
        (50, 32): ("red", "o", False),
    }
    points = []

    for min_bp, coverage in style_map:
        cohort_dir = root_dir / f"{min_bp}bp" / f"{coverage}x"
        pattern = f"{min_bp}bp_{coverage}x_*_alignment_distances_hapestry_kanpig.csv"
        for csv_path in sorted(cohort_dir.glob(pattern)):
            a = load_matrix(str(csv_path), delimiter=',')
            print(f"Loaded {csv_path} with shape {a.shape}. First 3 rows:\n{a[:3]}")
            count_h = np.sum(a[:, 1] < a[:, 3])
            count_k = np.sum(a[:, 3] < a[:, 1])
            total = a.shape[0]
            percent_h = 100.0 * count_h / total
            percent_k = 100.0 * count_k / total
            color, marker, filled = style_map[(min_bp, coverage)]
            points.append((percent_h, percent_k, color, marker, filled))

    if not points:
        raise ValueError(f"No matching CSV files found under {root_dir}")

    return points


def get_detail_csv(input_path: Path, sample_id: str) -> Path:
    """Resolve the CSV used by panels 2-4."""
    detail_candidates = input_path / "50bp" / "32x" / f"50bp_32x_{sample_id}_alignment_distances_hapestry_kanpig.csv"
    return detail_candidates


def load_all_cohort_csvs(root_dir: Path) -> dict:
    """Load all 4 cohort CSVs into a dict keyed by (min_bp, coverage)."""
    cohort_data = {}
    style_map = {
        (10, 8): ("blue", "o", True),
        (10, 32): ("red", "o", True),
        (50, 8): ("blue", "o", False),
        (50, 32): ("red", "o", False),
    }
    
    for min_bp, coverage in style_map:
        cohort_dir = root_dir / f"{min_bp}bp" / f"{coverage}x"
        pattern = f"{min_bp}bp_{coverage}x_*_alignment_distances_hapestry_kanpig.csv"
        csv_paths = sorted(cohort_dir.glob(pattern))
        if csv_paths:
            a = load_matrix(str(csv_paths[0]), delimiter=',')
            cohort_data[(min_bp, coverage)] = a
    
    return cohort_data


def log_bucket_medians(
    lengths: np.ndarray,
    values_h: np.ndarray,
    values_k: np.ndarray,
    n_bins: int = 12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-bucket medians for two value vectors using log-spaced length buckets."""
    valid = lengths > 0
    lengths = lengths[valid]
    values_h = values_h[valid]
    values_k = values_k[valid]

    if lengths.size == 0:
        raise ValueError("Column 0 must contain at least one positive value for log bucketing")

    min_len = float(np.min(lengths))
    max_len = float(np.max(lengths))
    if min_len == max_len:
        # Single-size dataset: return one bucket centered on that length.
        return np.array([min_len]), np.array([np.median(values_h)]), np.array([np.median(values_k)])

    edges = np.logspace(np.log10(min_len), np.log10(max_len), n_bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])

    bucket_idx = np.digitize(lengths, edges) - 1
    bucket_idx[bucket_idx == n_bins] = n_bins - 1

    med_h = []
    med_k = []
    kept_centers = []
    for i in range(n_bins):
        in_bucket = bucket_idx == i
        if np.any(in_bucket):
            kept_centers.append(centers[i])
            med_h.append(np.median(values_h[in_bucket]))
            med_k.append(np.median(values_k[in_bucket]))

    return np.array(kept_centers), np.array(med_h), np.array(med_k)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot distance-evaluation summaries from a CSV file or cohort root directory.")
    parser.add_argument("input_path", help="Path to an input CSV, or to a ROOT directory containing Xbp/Yx cohort CSVs")
    args = parser.parse_args()

    input_path = Path(args.input_path)
    left_panel_points = None
    if input_path.is_dir():
        left_panel_points = get_left_panel_points(input_path)

    detail_sample_id="HG002"
    detail_csv = get_detail_csv(input_path, detail_sample_id)
    a = load_matrix(str(detail_csv), delimiter=',')
    
    # Load all cohort data for multi-plot if directory input
    all_cohort_data = None
    if input_path.is_dir():
        all_cohort_data = load_all_cohort_csvs(input_path)

    fig, axes = plt.subplots(1, 2, figsize=(24, 5))
    grid_color = "0.9"
    right_grid_color = "0.9"

    # First plot: all samples.
    ax_left = axes[1]
    ax_left.set_aspect("equal", adjustable="box")
    ax_left.grid(True, color=grid_color, linewidth=0.8)
    if left_panel_points is None:
        count_h = np.sum(a[:, 1] < a[:, 2])
        count_k = np.sum(a[:, 2] < a[:, 1])
        total = a.shape[0]
        percent_h = 100.0 * count_h / total
        percent_k = 100.0 * count_k / total
        left_panel_points = [(percent_h, percent_k, "red", "o", False)]
    labels_used = set()
    style_labels = {
        ("blue", "o", True): "10bp 8x",
        ("red", "o", True): "10bp 32x",
        ("blue", "o", False): "50bp 8x",
        ("red", "o", False): "50bp 32x",
    }
    for percent_h, percent_k, color, marker, filled in left_panel_points:
        style_key = (color, marker, filled)
        label = style_labels.get(style_key) if style_key not in labels_used else None
        labels_used.add(style_key)
        mfc = color if filled else 'none'
        ms = 3 if filled else 8
        ax_left.plot(percent_h, percent_k, linestyle='None', marker=marker, color=color,
                 markerfacecolor=mfc, markersize=ms, label=label)
    ax_left.plot([0, 100], [0, 100], color=grid_color)
    ax_left.set_xlim(0, 30)
    ax_left.set_ylim(0, 30)
    ax_left.set_xlabel(r"% windows with $w_H < w_{TK}$")
    ax_left.set_ylabel(r"% windows with $w_{TK} < w_H$")
    ax_left.set_title("Fraction of windows in each sample")
    ax_left.legend(loc="upper left")

    # Second plot: 2D heatmap of point density using rows of A as (x, y).
    ax_right = axes[0]
    x = a[:, 1]
    y = a[:, 3]
    N_BINS = 100
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    xedges = np.linspace(x_min, x_max, N_BINS + 1)
    yedges = np.linspace(y_min, y_max, N_BINS + 1)
    counts, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    masked_counts = np.ma.masked_where(counts == 0, counts)
    im = ax_right.pcolormesh(xedges, yedges, masked_counts.T, cmap='viridis', shading='auto', zorder=2)
    #cbar = fig.colorbar(im, ax=ax_right)
    #cbar.set_label("Number of windows")

    low = min(x_min, y_min)
    high = max(x_max, y_max)
    ax_right.plot([low, high], [low, high], color=grid_color, linewidth=1, zorder=3)

    ax_right.set_aspect("equal", adjustable="box")
    ax_right.set_axisbelow(True)
    ax_right.grid(True, which="both", color=right_grid_color, linewidth=0.8)
    ax_right.set_xlabel("$w_H$")
    ax_right.set_ylabel("$w_{TK}$")
    ax_right.set_title(rf"Matching weight per window in {detail_sample_id} ($\geq$50bp, 32x)")

    # # Third panel: per-row difference (kanpig score - hapestry score).
    # ax_third = axes[2]
    # if all_cohort_data:
    #     cohort_order = [(10, 8), (10, 32), (50, 8), (50, 32)]
    #     cohort_labels = ["10bp 8x", "10bp 32x", "50bp 8x", "50bp 32x"]
    #     all_diffs = []
        
    #     for min_bp, coverage in cohort_order:
    #         if (min_bp, coverage) in all_cohort_data:
    #             cohort_a = all_cohort_data[(min_bp, coverage)]
    #             diff = cohort_a[:, 3] - cohort_a[:, 1]
    #             all_diffs.append(diff)
    #         else:
    #             all_diffs.append(np.array([]))
        
    #     rng = np.random.default_rng(0)
    #     for i, diff in enumerate(all_diffs):
    #         if len(diff) > 0:
    #             jitter = rng.uniform(-0.15, 0.15, size=len(diff))
    #             ax_third.scatter(i + jitter, diff, s=8, color='#4878D0', alpha=0.5, linewidths=0, zorder=2)
        
    #     ax_third.boxplot(all_diffs, positions=range(len(cohort_order)), widths=0.4, showfliers=False,
    #                      medianprops=dict(color='black', linewidth=1.5),
    #                      boxprops=dict(linewidth=1.2),
    #                      whiskerprops=dict(linewidth=1.2),
    #                      capprops=dict(linewidth=1.2),
    #                      zorder=3)
    #     ax_third.axhline(0, color=grid_color, linewidth=1)
    #     ax_third.set_xticks(range(len(cohort_order)))
    #     ax_third.set_xticklabels(cohort_labels)
    #     ax_third.set_ylabel(r"$w_{TK} - w_H$")
    #     ax_third.set_title(rf"Matching weight difference per window in {detail_sample_id} ($\geq$50bp, 32x)")
    # else:
    #     diff = a[:, 3] - a[:, 1]
    #     rng = np.random.default_rng(0)
    #     jitter = rng.uniform(-0.15, 0.15, size=len(diff))
    #     ax_third.scatter(jitter, diff, s=8, color='#4878D0', alpha=0.5, linewidths=0, zorder=2)
    #     ax_third.boxplot(diff, positions=[0], widths=0.3, showfliers=False,
    #                      medianprops=dict(color='black', linewidth=1.5),
    #                      boxprops=dict(linewidth=1.2),
    #                      whiskerprops=dict(linewidth=1.2),
    #                      capprops=dict(linewidth=1.2),
    #                      zorder=3)
    #     ax_third.axhline(0, color=grid_color, linewidth=1)
    #     ax_third.set_xticks([0])
    #     ax_third.set_xticklabels([rf"$\geq$50bp, 32x"])
    #     ax_third.set_ylabel(r"$w_{TK} - w_H$")
    #     ax_third.set_title(rf"Match weight difference per window in {detail_sample_id} ($\geq$50bp, 32x)")
    
    # ax_third.grid(True, axis='y', color=right_grid_color, linewidth=0.8)
    

    # # Fourth panel: median w_H and w_TK after log-bucketing by window length (column 0).
    # ax_fourth = axes[3]
    # lengths = a[:, 0]
    # centers, med_h, med_tk = log_bucket_medians(lengths, a[:, 1], a[:, 3], n_bins=12)

    # for x_i, y_tk, y_h in zip(centers, med_tk, med_h):
    #     ax_fourth.annotate(
    #         "",
    #         xy=(x_i, y_h),
    #         xytext=(x_i, y_tk),
    #         arrowprops=dict(arrowstyle="->", color="0.35", linewidth=1.0, alpha=0.9),
    #         zorder=2,
    #     )

    # ax_fourth.plot(centers, med_h, "o", color="#1b9e77", markersize=4, label=r"median($w_H$)")
    # ax_fourth.plot(centers, med_tk, "o", color="#d95f02", markersize=4, label=r"median($w_{TK}$)")
    # ax_fourth.set_xscale("log")
    # ax_fourth.set_xlabel("Window length bin (bp)")
    # ax_fourth.set_ylabel("Matching weight")
    # ax_fourth.set_title(rf"Matching weight median by length bucket in {detail_sample_id} ($\geq$50bp, 32x)")
    # ax_fourth.grid(True, which="major", color=right_grid_color, linewidth=0.8)
    # ax_fourth.legend(loc="best")

    fig.tight_layout(pad=0.15, w_pad=0.3)
    fig.subplots_adjust(left=0.04, right=0.995, bottom=0.11, top=0.93, wspace=0.12)
    plt.show()


if __name__ == "__main__":
    main()
