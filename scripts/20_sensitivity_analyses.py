"""
Re-implementation of four analyses that were reported in the manuscript but had
no retained script: SA1 (threshold sensitivity), SA6 (depth-matched diversity),
the post-hoc non-synonymous sensitivity analysis, and the mixed-infection screen.

IMPORTANT PROVENANCE NOTE
-------------------------
These four analyses were originally run ad hoc and the code was not retained, so
they could not be regenerated from the deposited repository. They are therefore
RE-IMPLEMENTED here from the pre-registration (Sections 9.5-9.6) and the methods
description, not recovered. Two consequences follow, and both are disclosed in the
manuscript and response letter:

  1. Values will differ from the originally published ones because the corrected
     strand-bias filter reclassified 16,067 of 289,054 concordant calls.
  2. Values may additionally differ for implementation reasons that cannot be
     excluded, because the original code no longer exists to diff against.

From this point the analyses are reproducible: this script is deposited.

Site definitions and codon classification are imported from the primary pipeline
so each analysis differs from the primary in exactly one respect.

Outputs:
  results/sa1_threshold_sensitivity.csv
  results/sa6_depth_matched.json
  results/sensitivity_all_nonsyn_detections.csv
  results/sensitivity_all_nonsyn_summary.csv
  results/mixed_infection_flagged.csv
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _pipeline import (  # noqa: E402
    SITE_DEFINITIONS, classify_variant, AF_THRESHOLD, DEPTH_THRESHOLD,
    parse_ivar, parse_lofreq, strand_bias_filter,
)

RESULTS_DIR = PROJECT_ROOT / "results"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "v1_20260525_prioritised.tsv"

THRESHOLD_CANDIDATES = [0.01, 0.02, 0.03, 0.05]
PRIMARY_AF = 0.03
CONSENSUS_AF = 0.50
MIXED_INFECTION_ISNV_CUTOFF = 50
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
SEED = 42

ADAPTIVE = ("adaptation", "reversion_to_wildtype")


def load():
    variants = pd.read_parquet(RESULTS_DIR / "corpus_variants.parquet")
    site_depths = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet").rename(
        columns={"run_accession": "sample"}
    )
    diversity = pd.read_parquet(RESULTS_DIR / "diversity_by_sample.parquet")
    segcov = pd.read_parquet(RESULTS_DIR / "segment_coverage_by_sample.parquet")
    return variants, site_depths, coverage, diversity, segcov


def _site_variant_rows(variants, site):
    """Concordant, strand-bias-passing calls at a site's three codon positions."""
    return variants[
        (variants["chrom"] == site["segment_acc"])
        & (variants["pos"].isin(site["nt_positions_1idx"]))
        & (variants["passes_strand_bias"])
    ]


# ---------------------------------------------------------------- SA1

def sa1_threshold_sensitivity(variants, site_depths):
    """Repeat the H1 prevalence analysis at each candidate AF threshold."""
    rows = []
    for threshold in THRESHOLD_CANDIDATES:
        for site in SITE_DEFINITIONS:
            sid = site["site_id"]
            adequate = set(
                site_depths[site_depths[f"{sid}_min_depth"] >= DEPTH_THRESHOLD]["sample"]
            )
            hits = _site_variant_rows(variants, site)
            hits = hits[(hits["af_mean"] >= threshold) & (hits["sample"].isin(adequate))]

            positive, subconsensus = set(), set()
            for v in hits.itertuples(index=False):
                idx = site["nt_positions_1idx"].index(v.pos)
                if classify_variant(v.ref, v.alt, idx, site["wt_codon"],
                                    site["wt_aa"], site["adapted_aa"],
                                    site["direction"]) in ADAPTIVE:
                    positive.add(v.sample)
                    if v.af_mean < CONSENSUS_AF:
                        subconsensus.add(v.sample)

            n = len(adequate)
            rows.append({
                "af_threshold": threshold,
                "site_id": sid,
                "gene": site["gene"],
                "aa_pos": site["aa_pos"],
                "denominator": n,
                "n_detected": len(positive),
                "n_subconsensus": len(subconsensus),
                "prevalence_pct": 100 * len(positive) / n if n else np.nan,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- SA6

def _nearest_neighbour_match(target_df, pool_df, feature_cols):
    """Match each target sample to its nearest unused pool sample by Euclidean depth profile."""
    pool = pool_df.reset_index(drop=True)
    pool_features = pool[feature_cols].to_numpy(dtype=float)
    used, matched = set(), []
    for row in target_df[feature_cols].to_numpy(dtype=float):
        distances = np.linalg.norm(pool_features - row, axis=1)
        order = np.argsort(distances)
        for idx in order:
            if idx not in used:
                used.add(int(idx))
                matched.append(pool.iloc[int(idx)])
                break
    return pd.DataFrame(matched)


def _compare(a_pi, b_pi):
    u, p = stats.mannwhitneyu(a_pi, b_pi, alternative="two-sided")
    r = 1 - (2 * u) / (len(a_pi) * len(b_pi))
    return float(u), float(p), float(r)


def sa6_depth_matched(diversity, segcov, coverage):
    """Repeat the diversity comparisons on depth-matched cattle subsets."""
    feature_cols = [f"{s}_covered" for s in SEGMENTS]
    df = (diversity.merge(segcov[["sample"] + feature_cols], on="sample")
                   .merge(coverage[["sample", "median_depth"]], on="sample"))

    out = {}
    cattle = df[df["host_category"] == "cattle"]

    for label, host, key in (("cattle_vs_feline_matched", "feline", "feline"),
                             ("cattle_vs_milk_matched", "retail_milk", "milk")):
        target = df[df["host_category"] == host]
        if target.empty:
            continue
        matched = _nearest_neighbour_match(target, cattle, feature_cols)
        u, p, r = _compare(matched["genome_pi"].to_numpy(),
                           target["genome_pi"].to_numpy())
        entry = {
            "n_cattle": int(len(matched)),
            f"n_{key}": int(len(target)),
            "u_statistic": u,
            "p_value": p,
            "rank_biserial_r": r,
        }
        if key == "feline":
            entry.update({
                "cattle_median_depth": float(matched["median_depth"].median()),
                "feline_median_depth": float(target["median_depth"].median()),
                "cattle_median_pi": float(matched["genome_pi"].median()),
                "feline_median_pi": float(target["genome_pi"].median()),
            })
        out[label] = entry
    return out


# ---------------------------------------------------------------- non-synonymous

def nonsyn_sensitivity(variants, site_depths, coverage):
    """Classify every variant at a Tier 1 codon, not just the pre-specified substitution."""
    host = dict(zip(coverage["sample"], coverage["host_category"]))
    detections, summary = [], []

    for site in SITE_DEFINITIONS:
        sid = site["site_id"]
        adequate = set(
            site_depths[site_depths[f"{sid}_min_depth"] >= DEPTH_THRESHOLD]["sample"]
        )
        hits = _site_variant_rows(variants, site)
        hits = hits[(hits["af_mean"] >= PRIMARY_AF) & (hits["sample"].isin(adequate))]

        ref_codon = site["wt_codon"]
        ref_aa = str(Seq(ref_codon).translate())
        n_pre = n_other = n_syn = 0
        positive_samples = set()

        for v in hits.itertuples(index=False):
            idx = site["nt_positions_1idx"].index(v.pos)
            mutant = list(ref_codon)
            mutant[idx] = v.alt
            alt_codon = "".join(mutant)
            alt_aa = str(Seq(alt_codon).translate())

            cls = classify_variant(v.ref, v.alt, idx, ref_codon, site["wt_aa"],
                                   site["adapted_aa"], site["direction"])
            if cls in ADAPTIVE:
                change_type = "pre-specified adaptation"
                n_pre += 1
            elif alt_aa != ref_aa:
                change_type = "other non-synonymous"
                n_other += 1
            else:
                change_type = "synonymous"
                n_syn += 1

            if change_type != "synonymous":
                positive_samples.add(v.sample)

            detections.append({
                "site_id": sid, "gene": site["gene"], "aa_position": site["aa_pos"],
                "ref_aa": ref_aa, "alt_aa": alt_aa,
                "aa_change": f"{ref_aa}{site['aa_pos']}{alt_aa}",
                "change_type": change_type, "sample": v.sample,
                "af_mean": v.af_mean, "is_subconsensus": bool(v.af_mean < CONSENSUS_AF),
                "host_category": host.get(v.sample, "unknown"),
                "nt_position": v.pos, "ref_nt": v.ref, "alt_nt": v.alt,
                "ref_codon": ref_codon, "alt_codon": alt_codon,
            })

        n = len(adequate)
        summary.append({
            "site_id": sid, "gene": site["gene"], "aa_position": site["aa_pos"],
            "denominator": n,
            "n_prespecified": n_pre, "n_other_nonsyn": n_other,
            "n_synonymous": n_syn, "n_any_nonsyn": n_pre + n_other,
            "pct_any_nonsyn": round(100 * (n_pre + n_other) / n, 3) if n else np.nan,
        })

    return pd.DataFrame(detections), pd.DataFrame(summary)


# ---------------------------------------------------------------- mixed infection

def mixed_infection_screen(variants, coverage):
    """Flag samples whose sub-consensus iSNV burden suggests a mixed infection."""
    manifest = pd.read_csv(MANIFEST, sep="\t").rename(columns={"run_accession": "sample"})
    bioproject = dict(zip(manifest["sample"], manifest["bioproject"]))

    sub = variants[
        (variants["passes_strand_bias"])
        & (variants["af_mean"] >= PRIMARY_AF)
        & (variants["af_mean"] < CONSENSUS_AF)
    ]
    counts = sub.groupby("sample").size().rename("n_isnv").reset_index()
    flagged = counts[counts["n_isnv"] > MIXED_INFECTION_ISNV_CUTOFF].copy()

    # Segment-distribution uniformity: normalised Shannon entropy over the 8 segments.
    acc_to_seg = {s["segment_acc"]: s["gene"] for s in SITE_DEFINITIONS}
    uniformity = {}
    for sample, grp in sub[sub["sample"].isin(flagged["sample"])].groupby("sample"):
        per_seg = grp.groupby("chrom").size()
        p = per_seg / per_seg.sum()
        h = -(p * np.log(p)).sum()
        uniformity[sample] = float(h / np.log(len(SEGMENTS))) if len(per_seg) > 1 else 0.0

    cov = coverage.set_index("sample")
    flagged["host_category"] = flagged["sample"].map(cov["host_category"])
    flagged["median_depth"] = flagged["sample"].map(cov["median_depth"])
    flagged["bioproject"] = flagged["sample"].map(bioproject)
    flagged["segment_uniformity"] = flagged["sample"].map(uniformity)
    return flagged.sort_values("sample").reset_index(drop=True)


def main():
    np.random.seed(SEED)
    variants, site_depths, coverage, diversity, segcov = load()

    print("SA1: threshold sensitivity...")
    sa1 = sa1_threshold_sensitivity(variants, site_depths)
    sa1.to_csv(RESULTS_DIR / "sa1_threshold_sensitivity.csv", index=False)
    pivot = sa1.pivot_table(index=["site_id", "gene", "aa_pos"],
                            columns="af_threshold", values="n_detected")
    print(pivot.to_string())

    print("\nSA6: depth-matched diversity...")
    sa6 = sa6_depth_matched(diversity, segcov, coverage)
    (RESULTS_DIR / "sa6_depth_matched.json").write_text(json.dumps(sa6, indent=2))
    print(json.dumps(sa6, indent=2))

    print("\nNon-synonymous sensitivity...")
    det, summ = nonsyn_sensitivity(variants, site_depths, coverage)
    det.to_csv(RESULTS_DIR / "sensitivity_all_nonsyn_detections.csv", index=False)
    summ.to_csv(RESULTS_DIR / "sensitivity_all_nonsyn_summary.csv", index=False)
    print(summ.to_string(index=False))

    print("\nMixed-infection screen...")
    flagged = mixed_infection_screen(variants, coverage)
    flagged.to_csv(RESULTS_DIR / "mixed_infection_flagged.csv", index=False)
    print(f"  flagged samples (>{MIXED_INFECTION_ISNV_CUTOFF} sub-consensus iSNVs): {len(flagged)}")
    if not flagged.empty:
        print(f"  hosts: {flagged['host_category'].value_counts().to_dict()}")
        print(f"  segment uniformity range: {flagged['segment_uniformity'].min():.2f}"
              f"–{flagged['segment_uniformity'].max():.2f}")


if __name__ == "__main__":
    main()
