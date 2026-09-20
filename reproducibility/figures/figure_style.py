"""Shared typography and export settings for ProvFold figures."""

from __future__ import annotations

import matplotlib as mpl

MM_PER_INCH = 25.4
DOUBLE_COLUMN_IN = 170 / MM_PER_INCH

BLUE = "#3B6FB6"
TEAL = "#278C82"
ORANGE = "#D97945"
PURPLE = "#8064A2"
DARK = "#202A35"
MID = "#69717A"
LIGHT = "#D8DEE5"
PALE = "#F3F5F7"


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 7.5,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 7.0,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            # Preserve the registered 170 mm canvas. Axes-level layout controls
            # provide the crop, avoiding a different physical width per figure.
            "savefig.bbox": None,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.13,
        1.03,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color=DARK,
    )


def clean_axes(ax, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=LIGHT, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
