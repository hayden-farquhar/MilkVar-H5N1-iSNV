"""
Figure 3: Corpus composition and QC summary (pre-registered Section 9.9).

Panel A: Sample counts by host category and BioProject.
Panel B: Per-sample median depth distribution by host category.
Panel C: Proportion of samples meeting depth threshold per Tier 1 site.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DEPTH_THRESHOLD = 100
HOST_ORDER = ["cattle", "feline", "retail_milk"]
HOST_DISPLAY = {"cattle": "Cattle", "feline": "Feline", "retail_milk": "Retail milk"}
HOST_COLOURS = {"cattle": "#4393C3", "feline": "#F4A582", "retail_milk": "#92C5DE"}

SITE_LABELS = [
    "PB2-627", "PB2-701", "PB2-591", "PB2-271", "PB2-631",
    "PA-497", "HA-226", "HA-228", "PB1F2-66", "NS1-92", "M2-31",
]
SITE_IDS = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11"]


def main():
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet")
    coverage = coverage.rename(columns={"run_accession": "sample"})

    manifest = pd.read_csv(
        DATA_DIR / "manifests" / "v1_20260525_all_hosts.tsv", sep="\t",
        usecols=["run_accession", "bioproject", "host_category"],
    )
    manifest = manifest.rename(columns={"run_accession": "sample"})
    manifest = manifest.drop_duplicates(subset="sample")

    site_depths = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")

    merged = coverage.merge(manifest[["sample", "bioproject"]], on="sample", how="left")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel A: Sample counts by host and BioProject
    ax = axes[0]
    bp_host = merged.groupby(["bioproject", "host_category"]).size().unstack(fill_value=0)
    bp_host = bp_host.reindex(columns=[h for h in HOST_ORDER if h in bp_host.columns], fill_value=0)
    top_bps = bp_host.sum(axis=1).nlargest(6).index
    bp_host_top = bp_host.loc[top_bps]

    bp_short = [bp.replace("PRJNA", "P") for bp in bp_host_top.index]
    x = np.arange(len(bp_short))
    width = 0.25
    for i, host in enumerate([h for h in HOST_ORDER if h in bp_host_top.columns]):
        ax.bar(x + i * width, bp_host_top[host], width,
               color=HOST_COLOURS[host], label=HOST_DISPLAY[host])
    ax.set_xticks(x + width)
    ax.set_xticklabels(bp_short, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel("Sample count", fontsize=10)
    ax.set_title("A. Samples by BioProject", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7)
    ax.set_yscale("log")

    # Panel B: Depth distribution by host
    ax = axes[1]
    data_b = []
    labels_b = []
    for host in HOST_ORDER:
        vals = coverage[coverage["host_category"] == host]["median_depth"].values
        if len(vals) > 0:
            data_b.append(vals)
            labels_b.append(HOST_DISPLAY[host])

    bp = ax.boxplot(data_b, tick_labels=labels_b, patch_artist=True, showfliers=False)
    for patch, host in zip(bp["boxes"], [h for h in HOST_ORDER if len(coverage[coverage["host_category"] == h]) > 0]):
        patch.set_facecolor(HOST_COLOURS[host])
        patch.set_alpha(0.7)

    for i, (host, vals) in enumerate(zip([h for h in HOST_ORDER if len(coverage[coverage["host_category"] == h]) > 0], data_b)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   alpha=0.15, s=5, color=HOST_COLOURS[host], zorder=3)

    ax.set_ylabel("Median depth (×)", fontsize=10)
    ax.set_title("B. Sequencing depth by host", fontsize=11, fontweight="bold")
    ax.axhline(y=DEPTH_THRESHOLD, color="red", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.text(0.98, DEPTH_THRESHOLD * 1.1, "100×", color="red", fontsize=7,
            ha="right", transform=ax.get_yaxis_transform())

    # Panel C: Coverage proportion per Tier 1 site
    ax = axes[2]
    proportions = []
    for sid in SITE_IDS:
        col = f"{sid}_min_depth"
        adequate = (site_depths[col] >= DEPTH_THRESHOLD).sum()
        proportions.append(adequate / len(site_depths))

    bars = ax.bar(range(len(SITE_LABELS)), proportions, color="#4393C3", alpha=0.8)
    ax.set_xticks(range(len(SITE_LABELS)))
    ax.set_xticklabels(SITE_LABELS, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("Proportion ≥ 100× depth", fontsize=10)
    ax.set_title("C. Coverage per Tier 1 site", fontsize=11, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.axhline(y=0.80, color="orange", linestyle="--", linewidth=0.8, alpha=0.7)

    for i, p in enumerate(proportions):
        ax.text(i, p + 0.02, f"{p:.0%}", fontsize=6, ha="center", va="bottom")

    plt.tight_layout()

    out_png = OUTPUT_DIR / "figure3_corpus_qc.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_png}")

    out_tiff = OUTPUT_DIR / "figure3_corpus_qc.tiff"
    fig.savefig(out_tiff, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_tiff}")

    plt.close(fig)


if __name__ == "__main__":
    main()
