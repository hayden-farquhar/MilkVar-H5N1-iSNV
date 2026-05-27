"""
Supplementary Figure S1: Genome-wide variant density profile by host category.

Sub-consensus iSNV density (50-nt windows, normalised per sample) across the
H5N1 genome, with Tier 1 site positions marked and segment boundaries indicated.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT_ROOT / "results"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SEGMENTS = {
    "PP755964.1": ("PB2", 2280, 0),
    "PP755963.1": ("PB1", 2274, 2280),
    "PP755962.1": ("PA", 2151, 4554),
    "PP755957.1": ("HA", 1704, 6705),
    "PP755960.1": ("NP", 1497, 8409),
    "PP755959.1": ("NA", 1410, 9906),
    "PP755958.1": ("M", 982, 11316),
    "PP755961.1": ("NS", 838, 12298),
}

HOST_COLOURS = {"cattle": "#4393C3", "feline": "#F4A582", "retail_milk": "#92C5DE"}
HOST_DISPLAY = {"cattle": "Cattle", "feline": "Feline", "retail_milk": "Retail milk"}
HOST_ORDER = ["cattle", "feline", "retail_milk"]

SITE_POSITIONS = {
    "PB2-627": 1879, "PB2-701": 2101, "PB2-591": 1771, "PB2-271": 811,
    "PB2-631": 1891, "PA-497": 4554 + 1489, "HA-226": 6705 + 712,
    "HA-228": 6705 + 718, "PB1F2-66": 2280 + 290, "NS1-92": 12298 + 274,
    "M2-31": 11316 + 779,
}

BIN_SIZE = 50


def main():
    variants = pd.read_parquet(RESULTS / "corpus_variants_subconsensus.parquet")
    cov = pd.read_parquet(RESULTS / "corpus_coverage.parquet")
    cov = cov.rename(columns={"run_accession": "sample"})
    variants = variants.merge(cov[["sample", "host_category"]], on="sample", how="left")

    variants = variants.copy()
    variants["genome_pos"] = variants.apply(
        lambda row: SEGMENTS[row["chrom"]][2] + row["pos"], axis=1
    )
    variants["bin"] = (variants["genome_pos"] // BIN_SIZE) * BIN_SIZE

    host_counts = cov["host_category"].value_counts().to_dict()

    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1, 0.5]})

    ax = axes[0]
    for host in HOST_ORDER:
        host_vars = variants[variants["host_category"] == host]
        bin_counts = host_vars.groupby("bin").size() / host_counts.get(host, 1)
        ax.fill_between(bin_counts.index, bin_counts.values,
                        alpha=0.4 if host == "cattle" else 0.6,
                        color=HOST_COLOURS[host], label=HOST_DISPLAY[host],
                        linewidth=0)
        ax.plot(bin_counts.index, bin_counts.values,
                color=HOST_COLOURS[host], linewidth=0.5, alpha=0.7)

    for acc, (gene, length, offset) in SEGMENTS.items():
        ax.axvline(x=offset, color="grey", linewidth=0.5, alpha=0.3)
        ax.text(offset + length / 2, ax.get_ylim()[1] * 0.95, gene,
                ha="center", va="top", fontsize=8, fontweight="bold", color="#333333")

    ax.set_ylabel("iSNVs per sample\n(50-nt windows)", fontsize=10)
    ax.set_title("Genome-wide sub-consensus iSNV density by host category",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")

    ax = axes[1]
    det = pd.read_csv(RESULTS / "adaptation_site_detections.csv")
    for site, gpos in SITE_POSITIONS.items():
        n_det = len(det[det["site_id"] == "S" + str(
            list(SITE_POSITIONS.keys()).index(site) + 1).zfill(2)])
        colour = "#DE2D26" if n_det > 0 else "#999999"
        ax.axvline(x=gpos, color=colour,
                   linewidth=1.5 if n_det > 0 else 0.8, alpha=0.7)
        ax.text(gpos, 0.5, site.split("-")[1] if "-" in site else site,
                rotation=90, fontsize=5, ha="center", va="center", color=colour)

    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Tier 1\nsites", fontsize=9)

    for acc, (gene, length, offset) in SEGMENTS.items():
        ax.axvline(x=offset, color="grey", linewidth=0.5, alpha=0.3)

    ax = axes[2]
    for acc, (gene, length, offset) in SEGMENTS.items():
        colour = "#4393C3" if list(SEGMENTS.keys()).index(acc) % 2 == 0 else "#92C5DE"
        ax.barh(0, length, left=offset, height=0.6, color=colour, alpha=0.5,
                edgecolor="white")
        ax.text(offset + length / 2, 0, gene, ha="center", va="center", fontsize=8)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Genome position (nt)", fontsize=10)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "figure_s1_variant_density.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_DIR / "figure_s1_variant_density.tiff", dpi=600,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved: outputs/figures/figure_s1_variant_density.{png,tiff}")


if __name__ == "__main__":
    main()
