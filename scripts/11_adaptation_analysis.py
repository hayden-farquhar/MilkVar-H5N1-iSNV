"""
Phase 5: Adaptation-site prevalence analysis (Section 9.4 of pre-registration).

For each of the 11 Tier 1 sites:
1. Compute site-level depth denominator (samples with >=100x at the specific position)
2. Identify samples with the adaptation-associated iSNV at AF >= 3% AND concordant
3. Compute prevalence with exact binomial CI and one-sided test
4. Stratify by host category
5. Report within-sample co-occurrence of multiple adaptation markers

Outputs:
  results/adaptation_site_prevalence.csv      — per-site prevalence table
  results/adaptation_site_detections.csv      — every detection with AF/depth
  results/adaptation_cooccurrence.csv         — multi-site positive samples
  results/adaptation_summary.json             — machine-readable summary
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from Bio import SeqIO
from Bio.Seq import Seq


PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data"

AF_THRESHOLD = 0.03
DEPTH_THRESHOLD = 100
BONFERRONI_K = 12
ALPHA = 0.05 / BONFERRONI_K  # 0.00417

SITE_DEFINITIONS = [
    {
        "site_id": "S01", "gene": "PB2", "aa_pos": 627,
        "segment_acc": "PP755964.1",
        "wt_aa": "E", "adapted_aa": "K",
        "wt_codon": "GAA", "direction": "forward",
        "nt_positions_1idx": [1879, 1880, 1881],
        "h3_label": None,
    },
    {
        "site_id": "S02", "gene": "PB2", "aa_pos": 701,
        "segment_acc": "PP755964.1",
        "wt_aa": "D", "adapted_aa": "N",
        "wt_codon": "GAC", "direction": "forward",
        "nt_positions_1idx": [2101, 2102, 2103],
        "h3_label": None,
    },
    {
        "site_id": "S03", "gene": "PB2", "aa_pos": 591,
        "segment_acc": "PP755964.1",
        "wt_aa": "Q", "adapted_aa": "K",
        "wt_codon": "CAA", "direction": "forward",
        "nt_positions_1idx": [1771, 1772, 1773],
        "h3_label": None,
    },
    {
        "site_id": "S04", "gene": "PB2", "aa_pos": 271,
        "segment_acc": "PP755964.1",
        "wt_aa": "T", "adapted_aa": "A",
        "wt_codon": "ACA", "direction": "forward",
        "nt_positions_1idx": [811, 812, 813],
        "h3_label": None,
    },
    {
        "site_id": "S05", "gene": "PB2", "aa_pos": 631,
        "segment_acc": "PP755964.1",
        "wt_aa": "M", "adapted_aa": "L",
        "wt_codon": "CTG", "direction": "reversion",
        "nt_positions_1idx": [1891, 1892, 1893],
        "h3_label": None,
    },
    {
        "site_id": "S06", "gene": "PA", "aa_pos": 497,
        "segment_acc": "PP755962.1",
        "wt_aa": "K", "adapted_aa": "R",
        "wt_codon": "AGG", "direction": "reversion",
        "nt_positions_1idx": [1489, 1490, 1491],
        "h3_label": None,
    },
    {
        "site_id": "S07", "gene": "HA", "aa_pos": 238,
        "segment_acc": "PP755957.1",
        "wt_aa": "Q", "adapted_aa": "L",
        "wt_codon": "CAA", "direction": "forward",
        "nt_positions_1idx": [712, 713, 714],
        "h3_label": "H3-226",
    },
    {
        "site_id": "S08", "gene": "HA", "aa_pos": 240,
        "segment_acc": "PP755957.1",
        "wt_aa": "G", "adapted_aa": "S",
        "wt_codon": "GGA", "direction": "forward",
        "nt_positions_1idx": [718, 719, 720],
        "h3_label": "H3-228",
    },
    {
        "site_id": "S09", "gene": "PB1-F2", "aa_pos": 66,
        "segment_acc": "PP755963.1",
        "wt_aa": "N", "adapted_aa": "S",
        "wt_codon": "AAT", "direction": "forward",
        "nt_positions_1idx": [290, 291, 292],
        "h3_label": None,
    },
    {
        "site_id": "S10", "gene": "NS1", "aa_pos": 92,
        "segment_acc": "PP755961.1",
        "wt_aa": "D", "adapted_aa": "E",
        "wt_codon": "GAC", "direction": "forward",
        "nt_positions_1idx": [274, 275, 276],
        "h3_label": None,
    },
    {
        "site_id": "S11", "gene": "M2", "aa_pos": 31,
        "segment_acc": "PP755958.1",
        "wt_aa": "S", "adapted_aa": "N",
        "wt_codon": "AGT", "direction": "forward",
        "nt_positions_1idx": [779, 780, 781],
        "h3_label": None,
    },
]


def load_data():
    variants = pd.read_parquet(RESULTS_DIR / "corpus_variants.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet")
    coverage = coverage.rename(columns={"run_accession": "sample"})
    site_depths = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")
    return variants, coverage, site_depths


def classify_variant(ref_nt: str, alt_nt: str, codon_pos_in_codon: int,
                     ref_codon: str, wt_aa: str, adapted_aa: str,
                     direction: str) -> str:
    """Determine if a variant at a codon position produces the adaptation change."""
    mutant_codon = list(ref_codon)
    mutant_codon[codon_pos_in_codon] = alt_nt
    mutant_codon = "".join(mutant_codon)
    try:
        mutant_aa = str(Seq(mutant_codon).translate())
    except Exception:
        return "unknown"

    if direction == "forward":
        if mutant_aa == adapted_aa:
            return "adaptation"
        elif mutant_aa != wt_aa:
            return "other_nonsynonymous"
        else:
            return "synonymous"
    else:  # reversion
        if mutant_aa == wt_aa:
            return "reversion_to_wildtype"
        elif mutant_aa != adapted_aa:
            return "other_nonsynonymous"
        else:
            return "synonymous"


def analyse_site(variants: pd.DataFrame, coverage: pd.DataFrame,
                 site_depths: pd.DataFrame, site: dict) -> tuple:
    """Run the full Section 9.4 analysis for one adaptation site."""
    acc = site["segment_acc"]
    positions = site["nt_positions_1idx"]
    ref_codon = site["wt_codon"] if site["direction"] == "forward" else None

    if site["direction"] == "reversion":
        ref_codon = site["wt_codon"]

    site_vars = variants[
        (variants["chrom"] == acc)
        & (variants["pos"].isin(positions))
    ].copy()

    depth_col = f"{site['site_id']}_min_depth"
    adequate = site_depths[site_depths[depth_col] >= DEPTH_THRESHOLD]
    denominator = len(adequate)
    adequate_samples = set(adequate["sample"])

    # Find adaptation-associated variants
    detections = []
    for _, var in site_vars.iterrows():
        if var["sample"] not in adequate_samples:
            continue
        if var["af_mean"] < AF_THRESHOLD:
            continue
        if not var["passes_strand_bias"]:
            continue

        # Determine which position in the codon this variant affects
        codon_pos = positions.index(var["pos"])

        # Get the current reference codon from the site definition
        current_ref_codon = site["wt_codon"]

        classification = classify_variant(
            var["ref"], var["alt"], codon_pos,
            current_ref_codon, site["wt_aa"], site["adapted_aa"],
            site["direction"]
        )

        if classification in ("adaptation", "reversion_to_wildtype"):
            # Sub-consensus detection (AF < 50%)
            is_subconsensus = var["af_mean"] < 0.50

            detections.append({
                "site_id": site["site_id"],
                "gene": site["gene"],
                "aa_position": site["aa_pos"],
                "wt_aa": site["wt_aa"],
                "adapted_aa": site["adapted_aa"],
                "direction": site["direction"],
                "classification": classification,
                "sample": var["sample"],
                "chrom": var["chrom"],
                "nt_position": var["pos"],
                "ref_nt": var["ref"],
                "alt_nt": var["alt"],
                "af_ivar": var["af_ivar"],
                "af_lofreq": var["af_lofreq"],
                "af_mean": var["af_mean"],
                "depth_ivar": var["depth_ivar"],
                "depth_lofreq": var["depth_lofreq"],
                "is_subconsensus": is_subconsensus,
            })

    detections_df = pd.DataFrame(detections)

    # For prevalence: count samples with at least one adaptation detection
    if not detections_df.empty:
        positive_samples = detections_df["sample"].unique()
        subconsensus_samples = detections_df[
            detections_df["is_subconsensus"]
        ]["sample"].unique()
    else:
        positive_samples = np.array([])
        subconsensus_samples = np.array([])

    numerator_all = len(positive_samples)
    numerator_sub = len(subconsensus_samples)

    # Exact binomial CI (Clopper-Pearson) and one-sided test
    def binomial_analysis(k, n):
        if n == 0:
            return 0.0, (0.0, 0.0), 1.0
        prevalence = k / n
        ci_low, ci_high = stats.binom.ppf(
            [ALPHA / 2, 1 - ALPHA / 2], n, prevalence
        ) / n if k > 0 else (0.0, 0.0)
        # Clopper-Pearson exact 95% CI.
        # The k == 0 and k == n boundary cases must use the same confidence level as
        # the interior case below (0.025 / 0.975). Using the Bonferroni-corrected
        # ALPHA here instead would silently mix two confidence levels within a single
        # column reported as a 95% CI, and would inflate the upper bound at exactly
        # the zero-count sites where that bound carries the most interpretive weight.
        if k == 0:
            ci_low = 0.0
            ci_high = 1 - 0.025 ** (1 / n)
        elif k == n:
            ci_low = 0.025 ** (1 / n)
            ci_high = 1.0
        else:
            ci_low = stats.beta.ppf(0.025, k, n - k + 1)
            ci_high = stats.beta.ppf(0.975, k + 1, n - k)
        # One-sided exact binomial test (H0: proportion = 0)
        p_value = 1.0 - stats.binom.cdf(k - 1, n, 1e-10) if k > 0 else 1.0
        # More precisely: P(X >= k) under H0: p -> 0, which is essentially 0 for k>=1
        # Use binomial test against a null of 0
        p_value = stats.binomtest(k, n, 1e-10, alternative="greater").pvalue if k > 0 else 1.0
        return prevalence, (ci_low, ci_high), p_value

    prev_all, ci_all, p_all = binomial_analysis(numerator_all, denominator)
    prev_sub, ci_sub, p_sub = binomial_analysis(numerator_sub, denominator)

    # Host-stratified results
    host_results = {}
    if not detections_df.empty:
        det_with_host = detections_df.merge(
            coverage[["sample", "host_category"]], on="sample", how="left"
        )
        for host in ["cattle", "feline", "retail_milk"]:
            host_adequate = adequate.merge(
                coverage[["sample", "host_category"]], on="sample", how="left"
            )
            host_denom = len(host_adequate[host_adequate["host_category"] == host])
            host_det = det_with_host[det_with_host["host_category"] == host]
            host_pos = host_det["sample"].nunique()
            host_sub = host_det[host_det["is_subconsensus"]]["sample"].nunique()

            if host_denom > 0:
                host_prev = host_pos / host_denom
                if host_pos > 0:
                    ci_l = stats.beta.ppf(0.025, host_pos, host_denom - host_pos + 1)
                    ci_h = stats.beta.ppf(0.975, host_pos + 1, host_denom - host_pos)
                else:
                    ci_l = 0.0
                    ci_h = 1 - 0.025 ** (1 / host_denom)
            else:
                host_prev = 0.0
                ci_l, ci_h = 0.0, 0.0

            host_results[host] = {
                "denominator": host_denom,
                "numerator_all": host_pos,
                "numerator_subconsensus": host_sub,
                "prevalence": host_prev,
                "ci_95": (ci_l, ci_h),
            }
    else:
        for host in ["cattle", "feline", "retail_milk"]:
            host_adequate = adequate.merge(
                coverage[["sample", "host_category"]], on="sample", how="left"
            )
            host_denom = len(host_adequate[host_adequate["host_category"] == host])
            host_results[host] = {
                "denominator": host_denom,
                "numerator_all": 0,
                "numerator_subconsensus": 0,
                "prevalence": 0.0,
                "ci_95": (0.0, 1 - 0.025 ** (1 / host_denom) if host_denom > 0 else 0.0),
            }

    h3_label = site.get("h3_label", "")
    site_label = f"{site['gene']}-{site['aa_pos']}"
    if h3_label:
        site_label += f" ({h3_label})"

    prevalence_row = {
        "site_id": site["site_id"],
        "gene": site["gene"],
        "aa_position": site["aa_pos"],
        "site_label": site_label,
        "wt_aa": site["wt_aa"],
        "adapted_aa": site["adapted_aa"],
        "direction": site["direction"],
        "denominator": denominator,
        "numerator_all": numerator_all,
        "numerator_subconsensus": numerator_sub,
        "prevalence_all": prev_all,
        "prevalence_subconsensus": prev_sub,
        "ci_95_low": ci_all[0],
        "ci_95_high": ci_all[1],
        "ci_sub_95_low": ci_sub[0],
        "ci_sub_95_high": ci_sub[1],
        "p_value": p_all,
        "significant_bonferroni": p_all < ALPHA,
        "max_af": detections_df["af_mean"].max() if not detections_df.empty else 0.0,
        "median_af": detections_df["af_mean"].median() if not detections_df.empty else 0.0,
    }

    for host, hr in host_results.items():
        prevalence_row[f"{host}_denom"] = hr["denominator"]
        prevalence_row[f"{host}_n_all"] = hr["numerator_all"]
        prevalence_row[f"{host}_n_sub"] = hr["numerator_subconsensus"]
        prevalence_row[f"{host}_prev"] = hr["prevalence"]

    return prevalence_row, detections_df


def analyse_cooccurrence(detections_df: pd.DataFrame) -> pd.DataFrame:
    """Identify samples positive at multiple Tier 1 sites."""
    if detections_df.empty:
        return pd.DataFrame()

    sites_per_sample = (
        detections_df.groupby("sample")["site_id"]
        .apply(lambda x: sorted(x.unique().tolist()))
        .reset_index()
    )
    sites_per_sample = sites_per_sample.copy()
    sites_per_sample["n_sites"] = sites_per_sample["site_id"].apply(len)
    multi = sites_per_sample[sites_per_sample["n_sites"] > 1].copy()
    multi.loc[:, "sites"] = multi["site_id"].apply(lambda x: ", ".join(x))

    if multi.empty:
        return pd.DataFrame()

    # Add host info
    sample_host = detections_df[["sample"]].drop_duplicates()
    # Will be merged later with coverage data

    return multi[["sample", "n_sites", "sites"]].sort_values("n_sites", ascending=False)


def main():
    print("Loading data...")
    variants, coverage, site_depths = load_data()
    print(f"  Variants: {len(variants):,} rows, {variants['sample'].nunique():,} samples")
    print(f"  Coverage: {len(coverage):,} samples")
    print(f"  Site depths: {len(site_depths):,} samples")
    print(f"  AF threshold: {AF_THRESHOLD*100:.0f}%")
    print(f"  Depth threshold: {DEPTH_THRESHOLD}x")
    print(f"  Bonferroni alpha: {ALPHA:.5f} (K={BONFERRONI_K})")
    print()

    all_prevalence = []
    all_detections = []

    for site in SITE_DEFINITIONS:
        print(f"Analysing {site['site_id']} {site['gene']}-{site['aa_pos']} "
              f"({site['wt_aa']}->{site['adapted_aa']}, {site['direction']})...")

        prev_row, det_df = analyse_site(variants, coverage, site_depths, site)
        all_prevalence.append(prev_row)
        if not det_df.empty:
            all_detections.append(det_df)

        n_sub = prev_row["numerator_subconsensus"]
        n_all = prev_row["numerator_all"]
        denom = prev_row["denominator"]
        print(f"  Denominator: {denom:,} samples with >=100x depth")
        print(f"  Detections (all AF): {n_all} ({prev_row['prevalence_all']*100:.2f}%)")
        print(f"  Detections (sub-consensus): {n_sub} ({prev_row['prevalence_subconsensus']*100:.2f}%)")
        if n_all > 0:
            print(f"  Max AF: {prev_row['max_af']:.3f}, Median AF: {prev_row['median_af']:.3f}")
            sig = "YES" if prev_row["significant_bonferroni"] else "no"
            print(f"  Significant (Bonferroni): {sig}")
        print()

    # Compile results
    prevalence_df = pd.DataFrame(all_prevalence)
    prevalence_df.to_csv(RESULTS_DIR / "adaptation_site_prevalence.csv", index=False)

    if all_detections:
        detections_df = pd.concat(all_detections, ignore_index=True)
        detections_df.to_csv(RESULTS_DIR / "adaptation_site_detections.csv", index=False)
    else:
        detections_df = pd.DataFrame()

    # Co-occurrence analysis
    print("=" * 60)
    print("CO-OCCURRENCE ANALYSIS")
    print("=" * 60)
    cooccurrence = analyse_cooccurrence(detections_df)
    if not cooccurrence.empty:
        cooccurrence = cooccurrence.merge(
            coverage[["sample", "host_category"]], on="sample", how="left"
        )
        cooccurrence.to_csv(RESULTS_DIR / "adaptation_cooccurrence.csv", index=False)
        print(f"  Samples positive at multiple sites: {len(cooccurrence)}")
        print(f"  Max sites in one sample: {cooccurrence['n_sites'].max()}")
        print()
        for _, row in cooccurrence.head(20).iterrows():
            print(f"  {row['sample']} ({row.get('host_category', '?')}): "
                  f"{row['n_sites']} sites — {row['sites']}")
    else:
        print("  No multi-site positive samples detected.")
    print()

    # Summary JSON
    summary = {
        "analysis_date": "2026-05-27",
        "af_threshold": AF_THRESHOLD,
        "depth_threshold": DEPTH_THRESHOLD,
        "bonferroni_k": BONFERRONI_K,
        "alpha_corrected": ALPHA,
        "total_samples": int(coverage["sample"].nunique()),
        "sites_analysed": len(SITE_DEFINITIONS),
        "sites_with_detections": int((prevalence_df["numerator_all"] > 0).sum()),
        "sites_significant": int(prevalence_df["significant_bonferroni"].sum()),
        "total_detections": len(detections_df),
        "total_subconsensus_detections": int(
            detections_df["is_subconsensus"].sum() if not detections_df.empty else 0
        ),
        "multi_site_samples": len(cooccurrence) if not cooccurrence.empty else 0,
    }
    with open(RESULTS_DIR / "adaptation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Sites with any detection: {summary['sites_with_detections']}/11")
    print(f"  Sites significant (Bonferroni): {summary['sites_significant']}/11")
    print(f"  Total adaptation detections: {summary['total_detections']}")
    print(f"  Sub-consensus detections: {summary['total_subconsensus_detections']}")
    print(f"  Multi-site positive samples: {summary['multi_site_samples']}")
    print()
    print("Output files:")
    print(f"  {RESULTS_DIR / 'adaptation_site_prevalence.csv'}")
    print(f"  {RESULTS_DIR / 'adaptation_site_detections.csv'}")
    print(f"  {RESULTS_DIR / 'adaptation_cooccurrence.csv'}")
    print(f"  {RESULTS_DIR / 'adaptation_summary.json'}")


if __name__ == "__main__":
    main()
