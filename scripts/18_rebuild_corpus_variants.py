"""
Rebuild results/corpus_variants.parquet with the corrected strand-bias filter.

Context: iVar reports REF_DP/ALT_DP as TOTAL depth per allele and REF_RV/ALT_RV as
the reverse-strand portion. Both merge_results.py and concordance_filter.py passed
the totals into the Fisher exact test as if they were forward-strand counts, which
inflated each forward cell by its own reverse count and left the strand-bias filter
too permissive. Both are now fixed; this script re-derives the concordant variant
table from the retained per-caller outputs so the corpus reflects the corrected filter.

Coverage statistics are NOT recomputed: corpus_coverage.parquet is derived from
samtools depth output and is unaffected by the strand-bias filter.

Outputs:
  results/corpus_variants.parquet          (overwritten)
  results/corpus_variants_3pct.parquet     (overwritten; AF >= 3% subset)
  results/corpus_variants_subconsensus.parquet (overwritten; 3% <= AF < 50%)
  results/strandfix_delta.json             (what changed vs the pre-fix table)
"""

import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from _pipeline import concordance  # noqa: E402  (uses the corrected strand-bias filter)

RESULTS_DIR = PROJECT_ROOT / "results"
VARIANTS_DIR = PROJECT_ROOT / "data" / "processed" / "vastai_results" / "variants"
BACKUP = RESULTS_DIR / "_pre_strandfix_backup" / "corpus_variants.parquet"

AF_THRESHOLD = 0.03
CONSENSUS_AF = 0.50


def main():
    samples = sorted(p.name[: -len(".ivar.tsv")] for p in VARIANTS_DIR.glob("*.ivar.tsv"))
    print(f"Rebuilding concordant variants for {len(samples)} samples...")

    try:
        from tqdm import tqdm
        iterator = tqdm(samples, unit="sample")
    except ImportError:
        iterator = samples

    frames = []
    for n, sample in enumerate(iterator, 1):
        df = concordance(sample, VARIANTS_DIR)
        if not df.empty:
            frames.append(df)
        if n % 1000 == 0 and frames:
            pd.concat(frames, ignore_index=True).to_parquet(
                RESULTS_DIR / "_corpus_variants_partial.parquet", index=False
            )

    variants = pd.concat(frames, ignore_index=True)
    variants.to_parquet(RESULTS_DIR / "corpus_variants.parquet", index=False)

    passing = variants[variants["passes_strand_bias"]]
    at3 = passing[passing["af_mean"] >= AF_THRESHOLD]
    at3.to_parquet(RESULTS_DIR / "corpus_variants_3pct.parquet", index=False)
    at3[at3["af_mean"] < CONSENSUS_AF].to_parquet(
        RESULTS_DIR / "corpus_variants_subconsensus.parquet", index=False
    )

    delta = {
        "n_concordant_calls": int(len(variants)),
        "n_pass_strand_bias_now": int(variants["passes_strand_bias"].sum()),
        "n_fail_strand_bias_now": int((~variants["passes_strand_bias"]).sum()),
        "n_at_3pct_now": int(len(at3)),
    }

    if BACKUP.exists():
        old = pd.read_parquet(BACKUP)
        key = ["sample", "chrom", "pos", "alt"]
        merged = old[key + ["passes_strand_bias"]].merge(
            variants[key + ["passes_strand_bias"]],
            on=key, how="outer", suffixes=("_old", "_new"), indicator=True,
        )
        both = merged[merged["_merge"] == "both"]
        newly_failing = int(
            (both["passes_strand_bias_old"] & ~both["passes_strand_bias_new"]).sum()
        )
        newly_passing = int(
            (~both["passes_strand_bias_old"] & both["passes_strand_bias_new"]).sum()
        )
        delta.update({
            "n_concordant_calls_before": int(len(old)),
            "n_pass_strand_bias_before": int(old["passes_strand_bias"].sum()),
            "verdict_flipped_to_fail": newly_failing,
            "verdict_flipped_to_pass": newly_passing,
            "pct_of_calls_reclassified": round(
                100 * (newly_failing + newly_passing) / max(len(both), 1), 4
            ),
        })

    (RESULTS_DIR / "strandfix_delta.json").write_text(json.dumps(delta, indent=2))
    partial = RESULTS_DIR / "_corpus_variants_partial.parquet"
    if partial.exists():
        partial.unlink()

    print(json.dumps(delta, indent=2))


if __name__ == "__main__":
    main()
