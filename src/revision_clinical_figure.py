"""Presentation-only clinical figure; reads saved estimates and performs no refit."""
from pathlib import Path
import argparse
import os

import numpy as np
import pandas as pd


def build_figure(source_root, revision_root, destination):
    source_root, revision_root, destination = Path(source_root), Path(revision_root), Path(destination)
    os.environ.setdefault("MPLCONFIGDIR", str(destination.parent / "matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    data = revision_root / "results"
    destination.parent.mkdir(parents=True, exist_ok=True)
    bands = pd.read_csv(data / "joint_auc_bands.csv")
    observations = pd.read_csv(source_root / "results" / "primary_intervals.csv", usecols=["auc_nmol_l"])
    assert len(observations) == int(bands.intervals.sum()) == 276
    curves = {name: pd.read_csv(data / f"fitted_curve_{name}.csv") for name in ("events", "unwarned")}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11.5,
                         "axes.labelsize": 11.5, "axes.titlesize": 12,
                         "xtick.labelsize": 11, "ytick.labelsize": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 5.1), sharex="col",
                             gridspec_kw={"height_ratios": [3.3, 1]})
    colors = ["#247BA0", "#B6465F"]
    low, high = observations.auc_nmol_l.min(), observations.auc_nmol_l.max()
    edges = np.linspace(low, high, 26)
    for col, (name, title, rate) in enumerate([
        ("events", "A  All events / never-warning", "total_per_100days"),
        ("unwarned", "B  Unwarned events / M15", "unwarned_per_100days")]):
        ax, hist = axes[0, col], axes[1, col]
        curve, color = curves[name], colors[col]
        assert np.all(curve.lower <= curve.rate_per_100days) and np.all(curve.upper >= curve.rate_per_100days)
        ax.plot(curve.auc_nmol_l, curve.rate_per_100days, color=color, linewidth=2)
        ax.fill_between(curve.auc_nmol_l, curve.lower, curve.upper, color=color, alpha=.17)
        ax.errorbar(bands.auc_median, bands[rate],
                    yerr=np.vstack([bands[rate]-bands[rate+"_lower"], bands[rate+"_upper"]-bands[rate]]),
                    fmt="o", color="#252525", markersize=5, linewidth=1.3, capsize=3.5)
        ax.set_title(title, loc="left", pad=8, fontweight="bold")
        ax.set_ylim(0, 16.6)
        ax.set_yticks([0, 4, 8, 12, 16])
        ax.grid(axis="y", alpha=.18, linewidth=.7)
        hist.hist(observations.auc_nmol_l, bins=edges, color=color, alpha=.75)
        for boundary in (.2, .4, .6):
            hist.axvline(boundary, color="#666666", linewidth=.8, linestyle="--")
        hist.set_xlim(low, high)
        hist.set_ylim(0, 38)
        hist.set_yticks([0, 20])
        hist.set_xticks([.5, 1, 1.5, 2])
        hist.tick_params(axis="x", pad=4)
    axes[0, 0].set_ylabel("Episodes per 100\nobserved-sensor days", labelpad=7)
    axes[1, 0].set_ylabel("Intervals", labelpad=11)
    fig.supxlabel("Earlier normalized MMTT C-peptide AUC (nmol/L)", fontsize=12, y=.025)
    handles = [Line2D([0], [0], color="#555555", linewidth=2, label="Adjusted mean"),
               Patch(facecolor="#888888", alpha=.20, edgecolor="none", label="95% cluster CI"),
               Line2D([0], [0], marker="o", color="#252525", linewidth=1.3, markersize=5, label="Band rate + 95% CI")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               fontsize=10.8, handlelength=2, columnspacing=1.4, bbox_to_anchor=(.54, 1.005))
    fig.subplots_adjust(left=.115, right=.985, top=.875, bottom=.14, wspace=.24, hspace=.12)
    fig.savefig(destination, dpi=240)
    plt.close(fig)
    print(str(destination))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root",type=Path,required=True)
    parser.add_argument("--revision-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    build_figure(args.source_root,args.revision_root,args.output)
