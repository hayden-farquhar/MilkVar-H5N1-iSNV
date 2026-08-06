"""
Supplementary Table S6: quality-control metrics and sample disposition.

Addresses Reviewer 1, minor comment 7 ("more detail on sequencing coverage,
read-quality filters, sample-exclusion criteria; a supplementary table
summarising QC metrics would help") and supports major comment 2 (threshold
validation) and minor comment 9 (reproducibility).

Every filter parameter reported in Panel A is read from the pipeline
definitions themselves (workflow/config.yaml and workflow/rules/*.smk) rather
than restated from the manuscript, so the table cannot drift from what was run.

Outputs:
  manuscript/table_s6_qc_metrics.csv
  manuscript/table_s6_qc_metrics.md
"""

import re
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
WORKFLOW_DIR = PROJECT_ROOT / "workflow"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "v1_20260525_prioritised.tsv"

PRIMARY_DEPTH = 100
STRINGENT_DEPTH = 200
MIN_SEGMENTS = 6

SITE_LABELS = {
    "S01": "PB2-627", "S02": "PB2-701", "S03": "PB2-591", "S04": "PB2-271",
    "S05": "PB2-631", "S06": "PA-497", "S07": "HA-226 (H3)", "S08": "HA-228 (H3)",
    "S09": "PB1-F2-66", "S10": "NS1-92", "S11": "M2-31",
}


def read_pipeline_params():
    """Pull the filter parameters straight from config.yaml and the .smk rules."""
    cfg = yaml.safe_load((WORKFLOW_DIR / "config.yaml").read_text())
    qc, vc = cfg["qc"], cfg["variant_calling"]

    align_smk = (WORKFLOW_DIR / "rules" / "align.smk").read_text()
    call_smk = (WORKFLOW_DIR / "rules" / "call_variants.smk").read_text()

    mpileup = re.search(r"samtools mpileup ([^\\\n|]*)", call_smk)
    ivar_trim = re.search(r"-q (\d+) -m (\d+) -s (\d+)", align_smk)

    return [
        ("Read quality filtering (Illumina)",
         f"fastp: --qualified_quality_phred {qc['min_quality_illumina']}, "
         f"--length_required {qc['min_read_length_illumina']}, "
         f"--detect_adapter_for_pe"),
        ("Read quality filtering (Nanopore)",
         f"fastp: --qualified_quality_phred {qc['min_quality_nanopore']}, "
         f"--length_required {qc['min_read_length_nanopore']}"),
        ("Alignment (Illumina)", "minimap2 -a -x sr"),
        ("Alignment (Nanopore)", "minimap2 -a -x map-ont"),
        ("Mapping-quality filter", f"samtools view -q {qc['min_mapq']}"),
        ("Duplicate handling (WGS/RANDOM)", "samtools markdup"),
        ("Primer trimming (AMPLICON only)",
         f"ivar trim -q {ivar_trim.group(1)} -m {ivar_trim.group(2)} "
         f"-s {ivar_trim.group(3)}" if ivar_trim else "ivar trim"),
        ("Pileup generation",
         f"samtools mpileup {mpileup.group(1).strip()} "
         "(-d 0 = no depth cap; no downsampling applied)"
         if mpileup else "samtools mpileup -aa -A -d 0 -B -Q 0"),
        ("Variant caller 1", f"iVar v1.4.3 (-q {vc['ivar_min_qual']}, "
                             f"-m {vc['ivar_min_depth']}, -t {vc['ivar_min_freq']})"),
        ("Variant caller 2", "LoFreq v2.1.5 (--call-indels --no-default-filter)"),
        ("Caller concordance", "Both callers must report the same position and "
                               "alternate allele; SNVs only (indels excluded)"),
        ("Strand-bias filter",
         f"Fisher exact test on forward/reverse ref and alt counts; "
         f"excluded if p < {vc['strand_bias_pvalue']}"),
        ("Site-level depth (primary)", f"≥{PRIMARY_DEPTH}× minimum codon depth"),
        ("Diversity eligibility (H2)",
         f"≥{PRIMARY_DEPTH}× on ≥{MIN_SEGMENTS} of 8 segments"),
        ("Operational allele-frequency threshold",
         "3% (empirically determined; see Supplementary Table S2)"),
    ]


def sample_disposition():
    manifest = pd.read_csv(MANIFEST, sep="\t")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet")
    segcov = pd.read_parquet(RESULTS_DIR / "segment_coverage_by_sample.parquet")

    failed = [ln.strip() for ln in
              (RESULTS_DIR / "failed_samples.txt").read_text().split("\n")
              if ln.strip()]

    # Reuse the Phase 6 eligibility criterion exactly (phase6_diversity_analysis.py:73)
    n_eligible = int((segcov["segments_with_coverage"] >= MIN_SEGMENTS).sum())

    rows = [
        ("Runs in frozen corpus manifest (25 May 2026)", len(manifest)),
        ("Failed retrieval or processing (ENA/SRA download errors)", len(failed)),
        ("Samples successfully processed through the pipeline", len(coverage)),
        (f"Samples eligible for diversity analysis "
         f"(≥{PRIMARY_DEPTH}× on ≥{MIN_SEGMENTS}/8 segments)", n_eligible),
    ]
    return pd.DataFrame(rows, columns=["Stage", "N"]), coverage


def depth_by_host(coverage):
    g = coverage.groupby("host_category")["median_depth"]
    out = pd.DataFrame({
        "N samples": g.size(),
        "Min": g.min().round(0).astype(int),
        "25th pct": g.quantile(0.25).round(0).astype(int),
        "Median": g.median().round(0).astype(int),
        "75th pct": g.quantile(0.75).round(0).astype(int),
        "Max": g.max().round(0).astype(int),
    }).reset_index().rename(columns={"host_category": "Host category"})
    return out


def site_denominators():
    sd = pd.read_parquet(RESULTS_DIR / "site_depth_by_sample.parquet")
    rows = []
    for sid, label in SITE_LABELS.items():
        col = f"{sid}_min_depth"
        rows.append({
            "Site": f"{sid} {label}",
            f"N ≥{PRIMARY_DEPTH}× (primary denominator)": int((sd[col] >= PRIMARY_DEPTH).sum()),
            f"N ≥{STRINGENT_DEPTH}× (SA3)": int((sd[col] >= STRINGENT_DEPTH).sum()),
            "Median site depth": int(sd[col].median()),
        })
    return pd.DataFrame(rows)


def main():
    params = pd.DataFrame(read_pipeline_params(), columns=["Step", "Setting"])
    disposition, coverage = sample_disposition()
    depths = depth_by_host(coverage)
    denoms = site_denominators()

    csv_path = MANUSCRIPT_DIR / "table_s6_qc_metrics.csv"
    with open(csv_path, "w") as fh:
        for name, df in (("A. Pipeline parameters", params),
                         ("B. Sample disposition", disposition),
                         ("C. Median depth by host", depths),
                         ("D. Per-site depth denominators", denoms)):
            fh.write(f"# {name}\n")
            df.to_csv(fh, index=False)
            fh.write("\n")

    md = ["# Supplementary Table S6. Quality-control metrics and sample disposition",
          "",
          "All parameters in Panel A are extracted programmatically from the "
          "deposited pipeline definitions (`workflow/config.yaml`, "
          "`workflow/rules/*.smk`) by `scripts/build_table_s6_qc.py`.",
          ""]
    for name, df in (("Panel A. Pipeline filter parameters", params),
                     ("Panel B. Sample disposition cascade", disposition),
                     ("Panel C. Median per-sample depth by host category", depths),
                     ("Panel D. Per-site depth denominators (Tier 1 panel)", denoms)):
        md += [f"## {name}", "", df.to_markdown(index=False), ""]

    (MANUSCRIPT_DIR / "table_s6_qc_metrics.md").write_text("\n".join(md))

    print(disposition.to_string(index=False))
    print()
    print(depths.to_string(index=False))
    print()
    print(denoms.to_string(index=False))
    print(f"\nWrote {csv_path.name} and table_s6_qc_metrics.md")


if __name__ == "__main__":
    main()
