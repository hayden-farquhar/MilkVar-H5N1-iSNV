"""
Figure 1: Adaptation-site heatmap (pre-registered Section 9.9).

Rows = samples (ordered by host category, then collection date).
Columns = 11 Tier 1 sites.
Cell colour = allele frequency of adaptation-associated variant (white → dark red).
Cells with site-level depth < 100x are greyed out.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DEPTH_THRESHOLD = 100

SITE_ORDER = [
    "S01", "S02", "S03", "S04", "S05",
    "S06", "S07", "S08", "S09", "S10", "S11",
]

SITE_LABELS = {
    "S01": "PB2\n627 E→K",
    "S02": "PB2\n701 D→N",
    "S03": "PB2\n591 Q→K",
    "S04": "PB2\n271 T→A",
    "S05": "PB2\n631 L→M*",
    "S06": "PA\n497 R→K*",
    "S07": "HA\n226 Q→L",
    "S08": "HA\n228 G→S",
    "S09": "PB1-F2\n66 N→S",
    "S10": "NS1\n92 D→E",
    "S11": "M2\n31 S→N",
}

HOST_ORDER = ["cattle", "feline", "retail_milk"]
HOST_COLOURS = {"cattle": "#4393C3", "feline": "#F4A582", "retail_milk": "#92C5DE"}
HOST_DISPLAY = {"cattle": "Cattle", "feline": "Feline", "retail_milk": "Retail milk"}


def load_data():
    detections = pd.read_csv(RESULTS_DIR / "adaptation_site_detections.csv")
    site_depths = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet")
    coverage = coverage.rename(columns={"run_accession": "sample"})

    manifest = pd.read_csv(
        DATA_DIR / "manifests" / "v1_20260525_all_hosts.tsv", sep="\t",
        usecols=["run_accession", "collection_date", "geo_loc_name"],
    )
    manifest = manifest.rename(columns={"run_accession": "sample"})
    manifest = manifest.drop_duplicates(subset="sample")

    return detections, site_depths, coverage, manifest


def build_matrix(detections, site_depths, coverage, manifest):
    samples_with_signal = set()
    for sid in SITE_ORDER:
        depth_col = f"{sid}_min_depth"
        adequate = site_depths[site_depths[depth_col] >= DEPTH_THRESHOLD]["sample"]
        samples_with_signal.update(adequate)

    det_positive = set(detections["sample"].unique())

    all_samples = site_depths.merge(
        coverage[["sample", "host_category"]], on="sample", how="left"
    ).merge(manifest[["sample", "collection_date"]], on="sample", how="left")

    positive_samples = all_samples[all_samples["sample"].isin(det_positive)].copy()

    positive_samples = positive_samples.copy()
    positive_samples["host_sort"] = positive_samples["host_category"].map(
        {h: i for i, h in enumerate(HOST_ORDER)}
    )
    positive_samples = positive_samples.sort_values(
        ["host_sort", "collection_date", "sample"]
    ).reset_index(drop=True)

    n = len(positive_samples)
    af_matrix = np.full((n, len(SITE_ORDER)), np.nan)
    depth_mask = np.zeros((n, len(SITE_ORDER)), dtype=bool)

    for j, sid in enumerate(SITE_ORDER):
        depth_col = f"{sid}_min_depth"
        for i, row in positive_samples.iterrows():
            if row[depth_col] >= DEPTH_THRESHOLD:
                depth_mask[i, j] = True

    site_det = detections.groupby(["sample", "site_id"])["af_mean"].max().reset_index()
    sample_to_idx = {s: i for i, s in enumerate(positive_samples["sample"])}

    for _, row in site_det.iterrows():
        if row["sample"] in sample_to_idx:
            i = sample_to_idx[row["sample"]]
            j = SITE_ORDER.index(row["site_id"])
            af_matrix[i, j] = row["af_mean"]

    af_matrix[~depth_mask] = np.nan
    af_matrix = np.where(depth_mask & np.isnan(af_matrix), 0.0, af_matrix)

    return positive_samples, af_matrix, depth_mask


def plot_heatmap(positive_samples, af_matrix, depth_mask):
    n_samples = len(positive_samples)
    n_sites = len(SITE_ORDER)

    fig_height = max(6, n_samples * 0.04 + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "adaptation", ["#FFFFFF", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"]
    )
    cmap.set_bad(color="#D9D9D9")

    masked_data = np.ma.masked_where(np.isnan(af_matrix), af_matrix)

    im = ax.pcolormesh(
        masked_data,
        cmap=cmap,
        vmin=0,
        vmax=max(0.5, np.nanmax(af_matrix) if not np.all(np.isnan(af_matrix)) else 0.5),
        edgecolors="white",
        linewidth=0.3,
    )

    ax.set_xticks(np.arange(n_sites) + 0.5)
    ax.set_xticklabels(
        [SITE_LABELS[s] for s in SITE_ORDER],
        fontsize=7, ha="center", va="top",
    )
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    host_boundaries = []
    prev_host = None
    for i, (_, row) in enumerate(positive_samples.iterrows()):
        if row["host_category"] != prev_host and prev_host is not None:
            host_boundaries.append(i)
        prev_host = row["host_category"]

    for boundary in host_boundaries:
        ax.axhline(y=boundary, color="black", linewidth=1.2)

    ax.set_yticks([])

    host_counts = positive_samples["host_category"].value_counts()
    host_midpoints = {}
    cum = 0
    for host in HOST_ORDER:
        count = host_counts.get(host, 0)
        if count > 0:
            host_midpoints[host] = cum + count / 2
            cum += count

    for host, mid in host_midpoints.items():
        ax.text(
            -0.3, mid, HOST_DISPLAY[host],
            ha="right", va="center", fontsize=9, fontweight="bold",
            color=HOST_COLOURS[host],
            transform=ax.get_yaxis_transform(),
        )

    cbar = fig.colorbar(im, ax=ax, shrink=0.5, aspect=20, pad=0.02)
    cbar.set_label("Allele frequency", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    legend_elements = [
        Patch(facecolor="#D9D9D9", edgecolor="grey", label="Depth < 100×"),
        Patch(facecolor="white", edgecolor="grey", label="AF = 0 (reference)"),
    ]
    ax.legend(
        handles=legend_elements, loc="lower right",
        fontsize=7, framealpha=0.9,
        bbox_to_anchor=(1.0, -0.08),
    )

    ax.set_xlabel("")
    ax.set_title(
        "Sub-consensus adaptation-associated iSNVs at Tier 1 sites",
        fontsize=11, fontweight="bold", pad=40,
    )

    note = "* Reversion sites (reference carries adapted allele; variant = reversion to wildtype)"
    fig.text(0.12, 0.01, note, fontsize=7, fontstyle="italic", color="#666666")

    plt.tight_layout()

    out_png = OUTPUT_DIR / "figure1_adaptation_heatmap.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_png}")

    out_tiff = OUTPUT_DIR / "figure1_adaptation_heatmap.tiff"
    fig.savefig(out_tiff, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_tiff}")

    plt.close(fig)


def main():
    print("Loading data...")
    detections, site_depths, coverage, manifest = load_data()
    print(f"  Detections: {len(detections)} rows across {detections['sample'].nunique()} samples")
    print(f"  Site depths: {len(site_depths)} samples")

    print("Building AF matrix...")
    positive_samples, af_matrix, depth_mask = build_matrix(
        detections, site_depths, coverage, manifest
    )
    print(f"  Samples with any detection: {len(positive_samples)}")
    print(f"  Matrix shape: {af_matrix.shape}")
    print(f"  Non-zero AF cells: {(af_matrix > 0).sum()}")
    print(f"  Insufficient depth cells (grey): {(~depth_mask).sum()}")

    print("Plotting heatmap...")
    plot_heatmap(positive_samples, af_matrix, depth_mask)
    print("Done.")


if __name__ == "__main__":
    main()
