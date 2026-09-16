#!/usr/bin/env python3
"""Generate a journal-style PRISMA flow figure from the registered counts."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


HERE = Path(__file__).resolve().parent
ROUND_DIR = HERE.parent / "ah_alpha_diversity"
SOURCE = ROUND_DIR / "systematic_review/prisma_counts.csv"
OUTPUT = Path(os.environ.get("PROVFOLD_FIGURE_OUTPUT", HERE / "generated_figures")).resolve()


def read_counts() -> dict[str, int]:
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        return {row["count_id"]: int(row["n"]) for row in csv.DictReader(handle)}


def box(ax, x: float, y: float, width: float, height: float, text: str) -> None:
    ax.add_patch(Rectangle((x, y), width, height, facecolor="white", edgecolor="#24364B", linewidth=0.8))
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=7.2, linespacing=1.25)


def arrow(ax, start, end) -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, linewidth=0.8, color="#24364B"))


def main() -> None:
    c = read_counts()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "figure.facecolor": "white", "axes.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(7.0, 7.3))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    box(ax, 0.08, 0.86, 0.52, 0.09, f"Records identified through databases\n(n = {c['identified']:,})")
    box(ax, 0.68, 0.86, 0.25, 0.09, f"Duplicates removed\n(n = {c['duplicates']:,})")
    arrow(ax, (0.60, 0.905), (0.68, 0.905))
    arrow(ax, (0.34, 0.855), (0.34, 0.785))
    box(ax, 0.08, 0.69, 0.52, 0.09, f"Titles and abstracts screened\n(n = {c['screened']:,})")
    box(ax, 0.68, 0.69, 0.25, 0.09, f"Records excluded\n(n = {c['screen_excluded']:,})")
    arrow(ax, (0.60, 0.735), (0.68, 0.735))
    arrow(ax, (0.34, 0.685), (0.34, 0.615))
    box(ax, 0.08, 0.52, 0.52, 0.09, f"Full-text reports sought\n(n = {c['sought']:,})")
    box(ax, 0.68, 0.52, 0.25, 0.09, f"Reports not retrieved\n(n = {c['not_retrieved']:,})")
    arrow(ax, (0.60, 0.565), (0.68, 0.565))
    arrow(ax, (0.34, 0.515), (0.34, 0.445))
    box(ax, 0.08, 0.34, 0.52, 0.10, f"Full-text reports assessed\n(n = {c['assessed']:,})")
    excluded = (
        f"Full-text reports excluded (n = {c['full_text_excluded']:,})\n"
        f"Wrong design: {c['wrong_design']}\nWrong population: {c['wrong_population']}\n"
        f"Wrong outcome or method: {c['wrong_outcome_method']}\n"
        f"Required data unavailable: {c['data_unavailable']}\nOther protocol incongruence: {c['protocol_incongruent']}"
    )
    box(ax, 0.65, 0.285, 0.32, 0.21, excluded)
    arrow(ax, (0.60, 0.39), (0.65, 0.39))
    arrow(ax, (0.34, 0.335), (0.34, 0.25))
    box(ax, 0.08, 0.13, 0.52, 0.115, f"Reports in qualitative AH review\n(n = {c['included']}; 24 comparison rows)\nReports with extractable AH–healthy alpha-diversity inputs: 6\nPrimary cohort families after deduplication: 5")
    ax.text(0.015, 0.90, "Identification", rotation=90, va="center", ha="center", fontweight="bold")
    ax.text(0.015, 0.70, "Screening", rotation=90, va="center", ha="center", fontweight="bold")
    ax.text(0.015, 0.42, "Eligibility", rotation=90, va="center", ha="center", fontweight="bold")
    ax.text(0.015, 0.19, "Included", rotation=90, va="center", ha="center", fontweight="bold")
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.04, top=0.99)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png", "svg", "tif"):
        fig.savefig(OUTPUT / f"Figure_S1.{suffix}", dpi=1000 if suffix == "tif" else (600 if suffix == "png" else None))
    plt.close(fig)


if __name__ == "__main__":
    main()
