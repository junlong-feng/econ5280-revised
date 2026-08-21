#!/usr/bin/env python3
"""Recreate the Lee (2008)-style electoral regression-discontinuity figure.

The course data contain the Democratic vote-share margin in election t (x)
and the Democratic vote share in election t+1 (y).  The original image is no
longer available, so this script creates a reproducible teaching illustration
using equal-width-bin means and separate fourth-order polynomial fits on the
two sides of the zero cutoff.
"""

import os
import tempfile
from pathlib import Path

# Set a writable cache before importing Matplotlib. This keeps the script
# portable in containers whose home directory is read-only.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "econ5280-matplotlib")
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "house.csv"
OUT_PDF = ROOT / "assets" / "Ch8RD_Lee_recreated.pdf"
OUT_PNG = ROOT / "assets" / "Ch8RD_Lee_recreated.png"


def binned_means(data: pd.DataFrame, width: float = 0.025) -> pd.DataFrame:
    """Return nonempty equal-width-bin means on [-1, 1]."""
    edges = np.arange(-1.0, 1.0 + width, width)
    work = data.copy()
    work["bin"] = pd.cut(
        work["x"], bins=edges, labels=False, include_lowest=True, right=False
    )
    means = work.groupby("bin", observed=True).agg(
        x=("x", "mean"), y=("y", "mean"), n=("y", "size")
    )
    return means.reset_index(drop=True)


def polynomial_curve(data: pd.DataFrame, side: str, degree: int = 4):
    """Fit and evaluate a polynomial separately to one side of the cutoff."""
    if side == "left":
        sample = data[data["x"] < 0]
        grid = np.linspace(-1.0, 0.0, 400, endpoint=False)
    else:
        sample = data[data["x"] >= 0]
        grid = np.linspace(0.0, 1.0, 400)
    coefficient = np.polyfit(sample["x"], sample["y"], deg=degree)
    return grid, np.polyval(coefficient, grid)


def main() -> None:
    data = pd.read_csv(DATA).dropna(subset=["x", "y"])
    bins = binned_means(data)
    left_x, left_y = polynomial_curve(data, "left")
    right_x, right_y = polynomial_curve(data, "right")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Palatino Linotype", "Palatino", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.scatter(
        bins["x"], bins["y"], s=20, facecolor="white", edgecolor="#2f4f6f",
        linewidth=0.8, zorder=3, label="Mean within vote-margin bin"
    )
    ax.plot(left_x, left_y, color="#a23b3b", linewidth=2.1)
    ax.plot(right_x, right_y, color="#a23b3b", linewidth=2.1,
            label="Separate fourth-order polynomial fits")
    ax.axvline(0, color="0.25", linestyle="--", linewidth=1.0)
    ax.annotate(
        "Democratic victory cutoff",
        xy=(0, 0.18), xytext=(0.08, 0.12),
        arrowprops={"arrowstyle": "->", "color": "0.25", "lw": 0.8},
        color="0.25", ha="left", va="center"
    )

    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Democratic vote-share margin in election $t$")
    ax.set_ylabel("Democratic vote share in election $t+1$")
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.grid(axis="y", color="0.90", linewidth=0.7)
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
