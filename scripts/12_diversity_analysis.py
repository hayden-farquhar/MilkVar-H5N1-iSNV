"""
Phase 6: Diversity comparisons (Section 9.5 of pre-registration).

1. Compute per-sample nucleotide diversity (pi), Shannon entropy, iSNV count.
2. Kruskal-Wallis test for genome-wide pi across host categories.
3. Pairwise Dunn's tests if omnibus significant.
4. Effect sizes (rank-biserial correlation) with BCa bootstrap CIs.
5. H3: Wilcoxon for iSNV count cattle vs retail-milk.
6. Generate Figures 2A, 2B, 2C.

Inputs:
  results/corpus_variants_3pct.parquet  — concordant variants at AF >= 3%
  results/segment_coverage_by_sample.parquet — per-segment position counts
  results/corpus_coverage.parquet — host category

Outputs:
  results/diversity_by_sample.parquet
  results/diversity_tests.json
  outputs/figures/figure2_diversity.png/.tiff
"""

import json
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AF_THRESHOLD = 0.03
DEPTH_THRESHOLD = 100
BONFERRONI_K = 12
ALPHA = 0.05 / BONFERRONI_K
MIN_SEGMENTS = 6
BOOTSTRAP_N = 10000
BOOTSTRAP_SEED = 42

SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

SEGMENT_ACCESSIONS = {
    "PP755964.1": "PB2",
    "PP755963.1": "PB1",
    "PP755962.1": "PA",
    "PP755957.1": "HA",
    "PP755960.1": "NP",
    "PP755959.1": "NA",
    "PP755958.1": "M",
    "PP755961.1": "NS",
}

HOST_ORDER = ["cattle", "feline", "retail_milk"]
HOST_DISPLAY = {"cattle": "Cattle", "feline": "Feline", "retail_milk": "Retail milk"}
HOST_COLOURS = {"cattle": "#4393C3", "feline": "#F4A582", "retail_milk": "#92C5DE"}


def load_data():
    variants = pd.read_parquet(RESULTS_DIR / "corpus_variants_3pct.parquet")
    seg_cov = pd.read_parquet(RESULTS_DIR / "segment_coverage_by_sample.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet")
    coverage = coverage.rename(columns={"run_accession": "sample"})
    return variants, seg_cov, coverage


def compute_diversity(variants: pd.DataFrame, seg_cov: pd.DataFrame,
                      coverage: pd.DataFrame) -> pd.DataFrame:
    """Compute pi, entropy, and iSNV count per sample."""
    eligible = seg_cov[seg_cov["segments_with_coverage"] >= MIN_SEGMENTS].copy()
    eligible_samples = set(eligible["sample"])

    sample_vars = variants[variants["sample"].isin(eligible_samples)].copy()
    sample_vars.loc[:, "segment"] = sample_vars["chrom"].map(SEGMENT_ACCESSIONS)

    results = []

    for sample in eligible["sample"]:
        sv = sample_vars[sample_vars["sample"] == sample]
        sc = eligible[eligible["sample"] == sample].iloc[0]

        pi_by_segment = {}
        entropy_values = []

        for gene in SEGMENT_ORDER:
            covered_positions = sc[f"{gene}_covered"]
            if covered_positions == 0:
                continue

            seg_vars = sv[sv["segment"] == gene]

            if seg_vars.empty:
                pi_by_segment[gene] = 0.0
                continue

            seg_pi_sum = 0.0
            for _, var in seg_vars.iterrows():
                af = var["af_mean"]
                site_pi = 2 * af * (1 - af)
                seg_pi_sum += site_pi

                if af > 0 and af < 1:
                    h = -(af * np.log(af) + (1 - af) * np.log(1 - af))
                    entropy_values.append(h)

            pi_by_segment[gene] = seg_pi_sum / covered_positions

        genome_pi_sum = 0.0
        for _, var in sv.iterrows():
            af = var["af_mean"]
            genome_pi_sum += 2 * af * (1 - af)

        genome_covered = sc["genome_covered"]
        genome_pi = genome_pi_sum / genome_covered if genome_covered > 0 else 0.0

        isnv_count = len(sv[(sv["af_mean"] >= AF_THRESHOLD) & (sv["af_mean"] < 0.50)])
        total_variant_count = len(sv)

        row = {
            "sample": sample,
            "genome_pi": genome_pi,
            "genome_covered_positions": genome_covered,
            "mean_entropy": np.mean(entropy_values) if entropy_values else 0.0,
            "isnv_count_subconsensus": isnv_count,
            "total_variant_count": total_variant_count,
        }

        for gene in SEGMENT_ORDER:
            row[f"{gene}_pi"] = pi_by_segment.get(gene, np.nan)

        results.append(row)

    diversity_df = pd.DataFrame(results)
    diversity_df = diversity_df.merge(
        coverage[["sample", "host_category"]], on="sample", how="left"
    )

    return diversity_df


def rank_biserial(x, y):
    """Rank-biserial correlation for Wilcoxon rank-sum test."""
    n1, n2 = len(x), len(y)
    u_stat = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return 1 - (2 * u_stat) / (n1 * n2)


def bootstrap_rank_biserial_ci(x, y, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    """BCa bootstrap CI for rank-biserial correlation."""
    rng = np.random.default_rng(seed)
    obs = rank_biserial(x, y)
    boot_vals = []
    for _ in range(n_boot):
        bx = rng.choice(x, size=len(x), replace=True)
        by = rng.choice(y, size=len(y), replace=True)
        try:
            boot_vals.append(rank_biserial(bx, by))
        except Exception:
            boot_vals.append(np.nan)
    boot_vals = np.array(boot_vals)
    boot_vals = boot_vals[~np.isnan(boot_vals)]

    # BCa correction
    z0 = stats.norm.ppf(np.mean(boot_vals < obs))

    jack_vals = []
    combined = np.concatenate([x, y])
    labels = np.concatenate([np.zeros(len(x)), np.ones(len(y))])
    for i in range(len(combined)):
        mask = np.ones(len(combined), dtype=bool)
        mask[i] = False
        jx = combined[mask & (labels == 0)]
        jy = combined[mask & (labels == 1)]
        if len(jx) > 0 and len(jy) > 0:
            jack_vals.append(rank_biserial(jx, jy))
    jack_vals = np.array(jack_vals)
    jack_mean = np.mean(jack_vals)
    a_hat = np.sum((jack_mean - jack_vals) ** 3) / (
        6 * (np.sum((jack_mean - jack_vals) ** 2) ** 1.5)
    ) if np.sum((jack_mean - jack_vals) ** 2) > 0 else 0

    alpha_levels = [0.025, 0.975]
    adjusted = []
    for al in alpha_levels:
        z_al = stats.norm.ppf(al)
        adj = stats.norm.cdf(z0 + (z0 + z_al) / (1 - a_hat * (z0 + z_al)))
        adjusted.append(adj)

    ci_low = np.nanpercentile(boot_vals, adjusted[0] * 100)
    ci_high = np.nanpercentile(boot_vals, adjusted[1] * 100)

    return obs, ci_low, ci_high


def run_statistical_tests(diversity_df: pd.DataFrame) -> dict:
    """Run all pre-registered statistical tests."""
    results = {}

    # H2: Kruskal-Wallis for genome-wide pi
    groups = {}
    for host in HOST_ORDER:
        vals = diversity_df[diversity_df["host_category"] == host]["genome_pi"].values
        if len(vals) > 0:
            groups[host] = vals

    group_values = [groups[h] for h in HOST_ORDER if h in groups]
    group_names = [h for h in HOST_ORDER if h in groups]

    if len(group_values) >= 2:
        kw_stat, kw_p = stats.kruskal(*group_values)
        results["h2_kruskal_wallis"] = {
            "statistic": float(kw_stat),
            "p_value": float(kw_p),
            "significant_bonferroni": kw_p < ALPHA,
            "groups": {h: {"n": len(groups[h]), "median": float(np.median(groups[h])),
                          "mean": float(np.mean(groups[h]))}
                      for h in group_names},
        }

        # Pairwise comparisons (always compute for reporting)
        pairwise = {}
        pairs = [("cattle", "feline"), ("cattle", "retail_milk"), ("feline", "retail_milk")]
        for h1, h2 in pairs:
            if h1 in groups and h2 in groups and len(groups[h1]) > 0 and len(groups[h2]) > 0:
                u_stat, u_p = stats.mannwhitneyu(
                    groups[h1], groups[h2], alternative="two-sided"
                )
                r_rb, ci_low, ci_high = bootstrap_rank_biserial_ci(
                    groups[h1], groups[h2]
                )
                pairwise[f"{h1}_vs_{h2}"] = {
                    "u_statistic": float(u_stat),
                    "p_value": float(u_p),
                    "rank_biserial_r": float(r_rb),
                    "r_ci_95": [float(ci_low), float(ci_high)],
                    "n1": len(groups[h1]),
                    "n2": len(groups[h2]),
                }
        results["h2_pairwise"] = pairwise

    # H3: Wilcoxon for iSNV count cattle vs retail-milk
    cattle_isnv = diversity_df[diversity_df["host_category"] == "cattle"]["isnv_count_subconsensus"].values
    milk_isnv = diversity_df[diversity_df["host_category"] == "retail_milk"]["isnv_count_subconsensus"].values

    if len(cattle_isnv) > 0 and len(milk_isnv) > 0:
        w_stat, w_p = stats.mannwhitneyu(
            cattle_isnv, milk_isnv, alternative="two-sided"
        )
        results["h3_wilcoxon_isnv"] = {
            "statistic": float(w_stat),
            "p_value": float(w_p),
            "significant_alpha05": w_p < 0.05,
            "cattle_median": float(np.median(cattle_isnv)),
            "milk_median": float(np.median(milk_isnv)),
            "cattle_n": len(cattle_isnv),
            "milk_n": len(milk_isnv),
        }

    # Shannon entropy comparison (same structure as pi)
    entropy_groups = {}
    for host in HOST_ORDER:
        vals = diversity_df[diversity_df["host_category"] == host]["mean_entropy"].values
        if len(vals) > 0:
            entropy_groups[host] = vals

    if len(entropy_groups) >= 2:
        ent_values = [entropy_groups[h] for h in HOST_ORDER if h in entropy_groups]
        ent_stat, ent_p = stats.kruskal(*ent_values)
        results["entropy_kruskal_wallis"] = {
            "statistic": float(ent_stat),
            "p_value": float(ent_p),
            "significant_bonferroni": ent_p < ALPHA,
            "groups": {h: {"n": len(entropy_groups[h]),
                          "median": float(np.median(entropy_groups[h]))}
                      for h in entropy_groups},
        }

    return results


def plot_figure2(diversity_df: pd.DataFrame):
    """Generate Figure 2: diversity comparison panels."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Genome-wide pi by host
    ax = axes[0]
    data_a = []
    positions_a = []
    for i, host in enumerate(HOST_ORDER):
        vals = diversity_df[diversity_df["host_category"] == host]["genome_pi"].values
        if len(vals) > 0:
            data_a.append(vals)
            positions_a.append(i)

    bp = ax.boxplot(data_a, positions=positions_a, widths=0.5,
                    patch_artist=True, showfliers=False)
    for i, (patch, host) in enumerate(zip(bp["boxes"], [h for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0])):
        patch.set_facecolor(HOST_COLOURS[host])
        patch.set_alpha(0.7)

    for i, (host, vals) in enumerate(zip([h for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0], data_a)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), positions_a[i]) + jitter, vals,
                   alpha=0.3, s=8, color=HOST_COLOURS[host], zorder=3)

    ax.set_xticks(positions_a)
    ax.set_xticklabels([HOST_DISPLAY[h] for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0],
                        fontsize=9)
    ax.set_ylabel("Nucleotide diversity (π)", fontsize=10)
    ax.set_title("A. Genome-wide π by host", fontsize=11, fontweight="bold")
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(-4, -4))

    # Panel B: Per-segment pi for PB2, PB1, PA
    ax = axes[1]
    long_segments = ["PB2", "PB1", "PA"]
    x_pos = 0
    tick_positions = []
    tick_labels = []
    for seg in long_segments:
        col = f"{seg}_pi"
        for i, host in enumerate(HOST_ORDER):
            vals = diversity_df[diversity_df["host_category"] == host][col].dropna().values
            if len(vals) > 0:
                bp = ax.boxplot([vals], positions=[x_pos], widths=0.6,
                                patch_artist=True, showfliers=False)
                bp["boxes"][0].set_facecolor(HOST_COLOURS[host])
                bp["boxes"][0].set_alpha(0.7)
            x_pos += 1
        tick_positions.append(x_pos - 2)
        tick_labels.append(seg)
        x_pos += 0.5

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_ylabel("Nucleotide diversity (π)", fontsize=10)
    ax.set_title("B. Per-segment π (longest segments)", fontsize=11, fontweight="bold")
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(-4, -4))

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=HOST_COLOURS[h], alpha=0.7,
                            label=HOST_DISPLAY[h]) for h in HOST_ORDER]
    ax.legend(handles=legend_elements, fontsize=7, loc="upper right")

    # Panel C: iSNV count by host
    ax = axes[2]
    data_c = []
    positions_c = []
    for i, host in enumerate(HOST_ORDER):
        vals = diversity_df[diversity_df["host_category"] == host]["isnv_count_subconsensus"].values
        if len(vals) > 0:
            data_c.append(vals)
            positions_c.append(i)

    bp = ax.boxplot(data_c, positions=positions_c, widths=0.5,
                    patch_artist=True, showfliers=False)
    for patch, host in zip(bp["boxes"], [h for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0]):
        patch.set_facecolor(HOST_COLOURS[host])
        patch.set_alpha(0.7)

    for i, (host, vals) in enumerate(zip([h for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0], data_c)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), positions_c[i]) + jitter, vals,
                   alpha=0.3, s=8, color=HOST_COLOURS[host], zorder=3)

    ax.set_xticks(positions_c)
    ax.set_xticklabels([HOST_DISPLAY[h] for h in HOST_ORDER if len(diversity_df[diversity_df["host_category"] == h]) > 0],
                        fontsize=9)
    ax.set_ylabel("Sub-consensus iSNV count", fontsize=10)
    ax.set_title("C. iSNV count by host", fontsize=11, fontweight="bold")

    plt.tight_layout()

    out_png = OUTPUT_DIR / "figure2_diversity.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_png}")

    out_tiff = OUTPUT_DIR / "figure2_diversity.tiff"
    fig.savefig(out_tiff, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_tiff}")

    plt.close(fig)


def main():
    print("Loading data...")
    variants, seg_cov, coverage = load_data()
    print(f"  Variants (3% AF): {len(variants):,} rows")
    print(f"  Segment coverage: {len(seg_cov):,} samples")

    eligible = seg_cov[seg_cov["segments_with_coverage"] >= MIN_SEGMENTS]
    print(f"  Eligible (>= {MIN_SEGMENTS} segments with coverage): {len(eligible):,}")
    print()

    print("Computing diversity metrics...")
    diversity_df = compute_diversity(variants, seg_cov, coverage)
    diversity_df.to_parquet(RESULTS_DIR / "diversity_by_sample.parquet", index=False)
    print(f"  Computed for {len(diversity_df)} samples")
    print()

    for host in HOST_ORDER:
        hd = diversity_df[diversity_df["host_category"] == host]
        if len(hd) > 0:
            print(f"  {HOST_DISPLAY[host]} (n={len(hd)}):")
            print(f"    Genome-wide pi: median={hd['genome_pi'].median():.2e}, "
                  f"mean={hd['genome_pi'].mean():.2e}")
            print(f"    Mean entropy: median={hd['mean_entropy'].median():.4f}")
            print(f"    Sub-consensus iSNVs: median={hd['isnv_count_subconsensus'].median():.0f}, "
                  f"mean={hd['isnv_count_subconsensus'].mean():.1f}")
    print()

    print("Running statistical tests...")
    test_results = run_statistical_tests(diversity_df)

    def convert_numpy(obj):
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    with open(RESULTS_DIR / "diversity_tests.json", "w") as f:
        json.dump(test_results, f, indent=2, default=convert_numpy)

    # Print key results
    if "h2_kruskal_wallis" in test_results:
        kw = test_results["h2_kruskal_wallis"]
        sig = "YES" if kw["significant_bonferroni"] else "no"
        print(f"  H2 Kruskal-Wallis: H={kw['statistic']:.2f}, p={kw['p_value']:.2e}, "
              f"significant (Bonferroni): {sig}")

    if "h2_pairwise" in test_results:
        for pair, res in test_results["h2_pairwise"].items():
            print(f"    {pair}: U={res['u_statistic']:.0f}, p={res['p_value']:.2e}, "
                  f"r_rb={res['rank_biserial_r']:.3f} "
                  f"[{res['r_ci_95'][0]:.3f}, {res['r_ci_95'][1]:.3f}]")

    if "h3_wilcoxon_isnv" in test_results:
        h3 = test_results["h3_wilcoxon_isnv"]
        sig = "YES" if h3["significant_alpha05"] else "no"
        print(f"  H3 Wilcoxon (cattle vs milk iSNV): W={h3['statistic']:.0f}, "
              f"p={h3['p_value']:.2e}, significant (alpha=0.05): {sig}")

    print()
    print("Generating Figure 2...")
    plot_figure2(diversity_df)
    print()
    print("Phase 6 complete.")
    print(f"Output files:")
    print(f"  {RESULTS_DIR / 'diversity_by_sample.parquet'}")
    print(f"  {RESULTS_DIR / 'diversity_tests.json'}")
    print(f"  {OUTPUT_DIR / 'figure2_diversity.png'}")
    print(f"  {OUTPUT_DIR / 'figure2_diversity.tiff'}")


if __name__ == "__main__":
    main()
