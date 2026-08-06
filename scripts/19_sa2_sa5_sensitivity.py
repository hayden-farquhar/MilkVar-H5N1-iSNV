"""
SA2-SA5: pre-registered sensitivity analyses (pre-registration Section 9.6, lines 420-426).

These four analyses were pre-registered and referenced in the manuscript Methods
("reported in Supplementary Table S5") but had never been computed. This script
produces them from retained artefacts, with no pipeline re-run:

  SA2  Single-caller analysis   -- iVar-only and LoFreq-only prevalence per site,
                                   plus the per-site concordance matrix
                                   (both / iVar only / LoFreq only).
                                   Requires the raw per-caller outputs, because
                                   corpus_variants.parquet retains concordant
                                   calls only.
  SA3  Stringent depth          -- repeat H1 at site depth >= 200x (primary: 100x).
  SA4  Platform stratification  -- repeat H1 split by Illumina / Oxford Nanopore.
  SA5  Library-type strat.      -- repeat H1 split by WGS / AMPLICON.

Site coordinates, codon classification, and the caller parsers / strand-bias
filter are imported from the primary pipeline so that each sensitivity analysis
differs from the primary in exactly one respect.

Outputs:
  results/sa2_concordance_matrix.csv   -- per-site caller concordance + prevalence
  results/sa2_sa5_sensitivity.csv      -- long-format table for Supplementary Table S5
  results/sa2_sa5_summary.json         -- machine-readable summary
  results/_cache_sa2_percaller.parquet -- intermediate per-caller calls at Tier 1 sites
"""

import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _pipeline import (  # noqa: E402
    SITE_DEFINITIONS, classify_variant, AF_THRESHOLD, DEPTH_THRESHOLD,
    parse_ivar, parse_lofreq, strand_bias_filter,
)

RESULTS_DIR = PROJECT_ROOT / "results"
VARIANTS_DIR = PROJECT_ROOT / "data" / "processed" / "vastai_results" / "variants"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "v1_20260525_prioritised.tsv"
CACHE = RESULTS_DIR / "_cache_sa2_percaller.parquet"

STRINGENT_DEPTH = 200
ADAPTIVE_CLASSES = ("adaptation", "reversion_to_wildtype")

# (chrom, pos) -> list of (site, codon_index)
SITE_POSITIONS = {}
for _site in SITE_DEFINITIONS:
    for _i, _p in enumerate(_site["nt_positions_1idx"]):
        SITE_POSITIONS.setdefault((_site["segment_acc"], _p), []).append((_site, _i))


def _is_adaptive(site, codon_idx, ref_nt, alt_nt):
    """True when this nucleotide change produces the site's adaptation/reversion AA."""
    return classify_variant(
        ref_nt, alt_nt, codon_idx,
        site["wt_codon"], site["wt_aa"], site["adapted_aa"], site["direction"],
    ) in ADAPTIVE_CLASSES


def _adaptive_hits(df, sample):
    """Filter a single caller's calls down to adaptive Tier 1 detections."""
    hits = []
    for row in df.itertuples(index=False):
        key = (row.chrom, row.pos)
        if key not in SITE_POSITIONS:
            continue
        if row.af < AF_THRESHOLD:
            continue
        if not strand_bias_filter({
            "ref_fwd": row.ref_fwd, "ref_rev": row.ref_rev,
            "alt_fwd": row.alt_fwd, "alt_rev": row.alt_rev,
        }):
            continue
        for site, codon_idx in SITE_POSITIONS[key]:
            if _is_adaptive(site, codon_idx, row.ref, row.alt):
                hits.append({
                    "sample": sample,
                    "site_id": site["site_id"],
                    "chrom": row.chrom,
                    "pos": row.pos,
                    "alt": row.alt,
                    "af": row.af,
                })
    return hits


def build_percaller_table(force=False):
    """Parse every raw per-caller output at Tier 1 positions. Cached (~9,100 files)."""
    if CACHE.exists() and not force:
        print(f"Using cached per-caller table: {CACHE.name}")
        return pd.read_parquet(CACHE)

    samples = sorted(p.name[: -len(".ivar.tsv")] for p in VARIANTS_DIR.glob("*.ivar.tsv"))
    print(f"Parsing per-caller outputs for {len(samples)} samples...")

    try:
        from tqdm import tqdm
        iterator = tqdm(samples, unit="sample")
    except ImportError:
        iterator = samples

    records = []
    for n, sample in enumerate(iterator, 1):
        ivar_path = VARIANTS_DIR / f"{sample}.ivar.tsv"
        lofreq_path = VARIANTS_DIR / f"{sample}.lofreq.vcf"
        if not lofreq_path.exists():
            continue
        try:
            ivar_df = parse_ivar(ivar_path)
            lofreq_df = parse_lofreq(lofreq_path)
        except Exception as exc:  # malformed/empty output for a failed run
            print(f"  WARNING: {sample} unparseable ({exc}); skipped")
            continue

        for caller, df in (("ivar", ivar_df), ("lofreq", lofreq_df)):
            if df.empty:
                continue
            for hit in _adaptive_hits(df, sample):
                hit["caller"] = caller
                records.append(hit)

        if n % 500 == 0:
            pd.DataFrame(records).to_parquet(CACHE, index=False)

    table = pd.DataFrame(records)
    table.to_parquet(CACHE, index=False)
    print(f"  {len(table)} adaptive Tier 1 calls across both callers -> {CACHE.name}")
    return table


def exact_test(k, n):
    """One-sided exact binomial test against H0: prevalence = 0, as in the primary."""
    if n == 0:
        return float("nan")
    return stats.binomtest(k, n, 1e-9, alternative="greater").pvalue


def run_sa2(percaller, site_depths, variants):
    """Caller-specific prevalence and the per-site concordance matrix.

    Two different notions of "concordant" appear here and must not be conflated:

      * n_primary — the primary analysis definition, reproduced exactly: mean AF
        across callers >= 3%, with the strand-bias test applied once to the iVar
        strand counts. This column reconciles with Table 2 by construction.
      * n_both — both callers independently call the variant at >= 3% of their OWN
        reported AF and each passes a strand-bias test on its OWN strand counts.
        This is stricter, and is the quantity the pre-registration's concordance
        matrix asks for.

    The two diverge mainly at very high depth, where LoFreq's unfiltered read counts
    give Fisher's exact test far more power than iVar's base-quality-filtered counts.
    """
    rows = []
    for site in SITE_DEFINITIONS:
        sid = site["site_id"]
        positions = site["nt_positions_1idx"]
        adequate = set(
            site_depths[site_depths[f"{sid}_min_depth"] >= DEPTH_THRESHOLD]["sample"]
        )
        n = len(adequate)

        sub = percaller[percaller["site_id"] == sid]
        ivar_pos = set(sub[sub["caller"] == "ivar"]["sample"]) & adequate
        lofreq_pos = set(sub[sub["caller"] == "lofreq"]["sample"]) & adequate

        both = ivar_pos & lofreq_pos
        ivar_only = ivar_pos - lofreq_pos
        lofreq_only = lofreq_pos - ivar_pos

        primary_hits = variants[
            (variants["chrom"] == site["segment_acc"])
            & (variants["pos"].isin(positions))
            & (variants["af_mean"] >= AF_THRESHOLD)
            & (variants["passes_strand_bias"])
            & (variants["sample"].isin(adequate))
        ]
        primary = {
            r.sample for r in primary_hits.itertuples(index=False)
            if _is_adaptive(site, positions.index(r.pos), r.ref, r.alt)
        }

        rows.append({
            "site_id": sid,
            "gene": site["gene"],
            "aa_position": site["aa_pos"],
            "denominator": n,
            "n_primary": len(primary),
            "n_both": len(both),
            "n_ivar_only": len(ivar_only),
            "n_lofreq_only": len(lofreq_only),
            "n_neither": n - len(ivar_pos | lofreq_pos),
            "prev_primary_pct": 100 * len(primary) / n if n else float("nan"),
            "prev_both_pct": 100 * len(both) / n if n else float("nan"),
            "prev_ivar_alone_pct": 100 * len(ivar_pos) / n if n else float("nan"),
            "prev_lofreq_alone_pct": 100 * len(lofreq_pos) / n if n else float("nan"),
            "p_primary": exact_test(len(primary), n),
        })
    return pd.DataFrame(rows)


def _prevalence_rows(variants, site_depths, analysis, stratum, samples, depth_min):
    """Concordant-call prevalence per site within one sample subset."""
    rows = []
    for site in SITE_DEFINITIONS:
        sid = site["site_id"]
        acc, positions = site["segment_acc"], site["nt_positions_1idx"]

        adequate = set(
            site_depths[site_depths[f"{sid}_min_depth"] >= depth_min]["sample"]
        ) & samples
        n = len(adequate)

        sv = variants[
            (variants["chrom"] == acc)
            & (variants["pos"].isin(positions))
            & (variants["af_mean"] >= AF_THRESHOLD)
            & (variants["passes_strand_bias"])
            & (variants["sample"].isin(adequate))
        ]

        positive = {
            r.sample for r in sv.itertuples(index=False)
            if _is_adaptive(site, positions.index(r.pos), r.ref, r.alt)
        }
        k = len(positive)
        rows.append({
            "analysis": analysis,
            "stratum": stratum,
            "site_id": sid,
            "gene": site["gene"],
            "aa_position": site["aa_pos"],
            "label": f"{site['gene']}-{site['aa_pos']} "
                     f"{site['wt_aa']}→{site['adapted_aa']}",
            "denominator": n,
            "detections": k,
            "prevalence_pct": 100 * k / n if n else float("nan"),
            "p_value": exact_test(k, n),
            "significant": (exact_test(k, n) < ALPHA) if n else False,
        })
    return rows


def main():
    site_depths = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")
    variants = pd.read_parquet(RESULTS_DIR / "corpus_variants.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet").rename(
        columns={"run_accession": "sample"}
    )
    manifest = pd.read_csv(MANIFEST, sep="\t").rename(
        columns={"run_accession": "sample"}
    )
    all_samples = set(site_depths["sample"])

    # --- SA2 -------------------------------------------------------------
    print("\n== SA2: single-caller analysis ==")
    percaller = build_percaller_table()
    sa2 = run_sa2(percaller, site_depths, variants)
    sa2.to_csv(RESULTS_DIR / "sa2_concordance_matrix.csv", index=False)
    print(sa2[["site_id", "gene", "aa_position", "n_primary", "n_both",
               "n_ivar_only", "n_lofreq_only"]].to_string(index=False))

    long_rows = []

    # --- SA3 -------------------------------------------------------------
    print(f"\n== SA3: stringent depth >= {STRINGENT_DEPTH}x ==")
    long_rows += _prevalence_rows(variants, site_depths, "SA3",
                                  f">={DEPTH_THRESHOLD}x (primary)",
                                  all_samples, DEPTH_THRESHOLD)
    long_rows += _prevalence_rows(variants, site_depths, "SA3",
                                  f">={STRINGENT_DEPTH}x (stringent)",
                                  all_samples, STRINGENT_DEPTH)

    # --- SA4 -------------------------------------------------------------
    print("== SA4: platform stratification ==")
    for platform, grp in coverage.groupby("platform"):
        long_rows += _prevalence_rows(variants, site_depths, "SA4", platform,
                                      set(grp["sample"]), DEPTH_THRESHOLD)

    # --- SA5 -------------------------------------------------------------
    print("== SA5: library-strategy stratification ==")
    for strategy, grp in manifest.groupby("library_strategy"):
        long_rows += _prevalence_rows(variants, site_depths, "SA5", strategy,
                                      set(grp["sample"]) & all_samples,
                                      DEPTH_THRESHOLD)

    long = pd.DataFrame(long_rows)
    long.to_csv(RESULTS_DIR / "sa2_sa5_sensitivity.csv", index=False)

    summary = {
        "af_threshold": AF_THRESHOLD,
        "alpha_bonferroni": ALPHA,
        "primary_depth": DEPTH_THRESHOLD,
        "stringent_depth": STRINGENT_DEPTH,
        "sa2": {
            "n_samples_parsed": int(percaller["sample"].nunique()),
            "total_ivar_only": int(sa2["n_ivar_only"].sum()),
            "total_lofreq_only": int(sa2["n_lofreq_only"].sum()),
            "total_both": int(sa2["n_both"].sum()),
        },
        "strata": {
            a: sorted(g["stratum"].unique().tolist())
            for a, g in long.groupby("analysis")
        },
    }
    (RESULTS_DIR / "sa2_sa5_summary.json").write_text(json.dumps(summary, indent=2))

    print("\nWrote:")
    for f in ("sa2_concordance_matrix.csv", "sa2_sa5_sensitivity.csv",
              "sa2_sa5_summary.json"):
        print(f"  results/{f}")


if __name__ == "__main__":
    main()
