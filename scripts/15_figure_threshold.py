"""
Figure 4: Threshold determination (pre-registered Section 9.9).

Caller concordance rate, strand-bias exclusion rate, and replicate concordance
rate as a function of AF threshold (1-5%), with the chosen operational
threshold (3%) highlighted.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

THRESHOLDS = [1, 2, 3, 5]

# C1 and C4 are properties of the Phase 4 validation harness (caller concordance
# and synthetic spike-in sensitivity), not of the corpus variant table.
CONCORDANCE = [45.5, 62.9, 72.2, 80.5]
SPIKE_SENSITIVITY = [93.3, 93.3, 93.3, 46.7]  # V3 95:5 mixture


def strand_bias_exclusion_rates():
    """C2, derived from the corpus variant table rather than transcribed.

    These values were previously hardcoded from PROGRESS.md and had drifted from
    the data they described. Computing them here keeps the figure consistent with
    the deposited variant table, including after the strand-bias filter correction.
    """
    import pandas as pd

    variants = pd.read_parquet(PROJECT_ROOT / "results" / "corpus_variants.parquet")
    return [
        round(100 * (~variants[variants["af_mean"] >= t / 100]["passes_strand_bias"]).mean(), 1)
        for t in THRESHOLDS
    ]


STRAND_BIAS = strand_bias_exclusion_rates()


def main():
    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(THRESHOLDS, CONCORDANCE, "o-", color="#4393C3", linewidth=2,
             markersize=8, label="Caller concordance (C1)", zorder=5)
    ax1.plot(THRESHOLDS, SPIKE_SENSITIVITY, "s-", color="#F4A582", linewidth=2,
             markersize=8, label="Spike-in sensitivity (C4)", zorder=5)
    ax1.plot(THRESHOLDS, STRAND_BIAS, "^-", color="#999999", linewidth=2,
             markersize=8, label="Strand-bias exclusion (C2)", zorder=5)

    ax1.axhline(y=80, color="#4393C3", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.text(5.1, 81, "C1 gate (80%)", fontsize=7, color="#4393C3", va="bottom")

    ax1.axhline(y=80, color="#F4A582", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.text(5.1, 78, "C4 gate (80%)", fontsize=7, color="#F4A582", va="top")

    ax1.axhline(y=20, color="#999999", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.text(5.1, 21, "C2 gate (20%)", fontsize=7, color="#999999", va="bottom")

    # Highlight chosen threshold
    ax1.axvline(x=3, color="#DE2D26", linestyle="-", linewidth=2, alpha=0.3, zorder=1)
    ax1.annotate("Chosen threshold: 3%\n(contingency activation)",
                 xy=(3, 72.2), xytext=(3.5, 55),
                 fontsize=9, fontweight="bold", color="#DE2D26",
                 arrowprops=dict(arrowstyle="->", color="#DE2D26", lw=1.5),
                 ha="left")

    ax1.set_xlabel("AF threshold (%)", fontsize=11)
    ax1.set_ylabel("Rate (%)", fontsize=11)
    ax1.set_title("Threshold determination: four-criterion evaluation",
                  fontsize=12, fontweight="bold")
    ax1.set_xticks(THRESHOLDS)
    ax1.set_xlim(0.5, 5.5)
    ax1.set_ylim(0, 105)
    ax1.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3,
               frameon=False)
    ax1.grid(True, alpha=0.2)

    plt.tight_layout()

    out_png = OUTPUT_DIR / "figure4_threshold.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_png}")

    out_tiff = OUTPUT_DIR / "figure4_threshold.tiff"
    fig.savefig(out_tiff, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_tiff}")

    plt.close(fig)


if __name__ == "__main__":
    main()
