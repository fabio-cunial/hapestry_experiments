#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def load_matrix(path: str, delimiter: str = ',') -> np.ndarray:
    """Load a numeric matrix from a text file with a specified delimiter."""
    data = np.loadtxt(path, delimiter=delimiter)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


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
    fig, axes = plt.subplots(1, 4, figsize=(24, 5))
    grid_color = "0.8"
    right_grid_color = "0.9"

    ax_left = axes[0]
    ax_left.set_aspect("equal", adjustable="box")
    ax_left.grid(True, color=grid_color, linewidth=0.8)

    a = load_matrix("join_cleaned.csv")
    count_h = np.sum(a[:, 1] < a[:, 2])
    count_k = np.sum(a[:, 2] < a[:, 1])
    total = a.shape[0]
    percent_h = 100.0 * count_h / total
    percent_k = 100.0 * count_k / total
    ax_left.plot(percent_h, percent_k, "ob", markerfacecolor='none', markersize=12, label=r"$\geq$50bp, 32x")
    ax_left.plot([0, 100], [0, 100], color=grid_color)
    ax_left.set_xlabel(r"% windows with $d_H < d_K$")
    ax_left.set_ylabel(r"% windows with $d_K < d_H$")
    ax_left.set_title("% windows in each sample")
    ax_left.legend(loc="upper left")

    # Right panel: 2D heatmap of point density using rows of A as (x, y).
    ax_right = axes[1]
    a = load_matrix("join_cleaned.csv", delimiter=',')
    x = a[:, 1]
    y = a[:, 2]
    N_BINS = 100
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    xedges = np.linspace(x_min, x_max, N_BINS + 1)
    yedges = np.linspace(y_min, y_max, N_BINS + 1)
    counts, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    masked_counts = np.ma.masked_where(counts == 0, counts)
    im = ax_right.pcolormesh(xedges, yedges, masked_counts.T, cmap='viridis', shading='auto', zorder=2)
    cbar = fig.colorbar(im, ax=ax_right)
    cbar.set_label("Number of windows")

    low = min(x_min, y_min)
    high = max(x_max, y_max)
    ax_right.plot([low, high], [low, high], color=grid_color, linewidth=1, zorder=3)

    ax_right.set_aspect("equal", adjustable="box")
    ax_right.set_axisbelow(True)
    ax_right.grid(True, which="both", color=right_grid_color, linewidth=0.8)
    ax_right.set_xlabel("$d_H$")
    ax_right.set_ylabel("$d_K$")
    ax_right.set_title(r"Pair edit distance of every window in HG002 ($\geq$50bp, 32x)")

    # Third panel: per-row difference (kanpig score - hapestry score).
    ax_third = axes[2]
    a = load_matrix("join_cleaned.csv", delimiter=',')
    diff = a[:, 2] - a[:, 1]
    rng = np.random.default_rng(0)
    jitter = rng.uniform(-0.15, 0.15, size=len(diff))
    ax_third.scatter(jitter, diff, s=8, color='#4878D0', alpha=0.5, linewidths=0, zorder=2)
    ax_third.boxplot(diff, positions=[0], widths=0.3, showfliers=False,
                     medianprops=dict(color='black', linewidth=1.5),
                     boxprops=dict(linewidth=1.2),
                     whiskerprops=dict(linewidth=1.2),
                     capprops=dict(linewidth=1.2),
                     zorder=3)
    ax_third.axhline(0, color=grid_color, linewidth=1)
    ax_third.set_xticks([0])
    ax_third.set_xticklabels([r"$\geq$50bp, 32x"])
    ax_third.set_ylabel(r"$d_K - d_H$")
    ax_third.set_title("Pair edit distance difference of every window in HG002")
    ax_third.grid(True, axis='y', color=right_grid_color, linewidth=0.8)

    # Fourth panel: median d_H and d_K after log-bucketing by window length (column 0).
    ax_fourth = axes[3]
    a = load_matrix("join_cleaned.csv", delimiter=',')
    lengths = a[:, 0]
    centers, med_h, med_k = log_bucket_medians(lengths, a[:, 1], a[:, 2], n_bins=12)

    for x_i, y_k, y_h in zip(centers, med_k, med_h):
        ax_fourth.annotate(
            "",
            xy=(x_i, y_h),
            xytext=(x_i, y_k),
            arrowprops=dict(arrowstyle="->", color="0.35", linewidth=1.0, alpha=0.9),
            zorder=2,
        )

    ax_fourth.plot(centers, med_h, "o", color="#1b9e77", markersize=4, label=r"median($d_H$)")
    ax_fourth.plot(centers, med_k, "o", color="#d95f02", markersize=4, label=r"median($d_K$)")
    ax_fourth.set_xscale("log")
    ax_fourth.set_xlabel("Window length bin (bp)")
    ax_fourth.set_ylabel("Pair edit distance")
    ax_fourth.set_title("Median by bucket (log length)")
    ax_fourth.grid(True, which="major", color=right_grid_color, linewidth=0.8)
    ax_fourth.legend(loc="best")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
