#!/usr/bin/env python3
"""
phase5b_genotype_stratification.py

Revision (R1) analysis addressing Virology reviewers' PB2-631 "reversion" concern.

Both reviewers noted that the D1.1 genotype never acquired the PB2 631L mammalian
adaptation, so a 631 L->M call in a D1.1-genotype sample aligned to the B3.13
reference (which carries 631L) is NOT a within-host reversion: it is the ancestral
avian state of a different genotype, i.e. a reference-relative genotype difference.

This script re-classifies every PB2-631 and PB2-701 detection by genotype
(BioProject proxy: PRJNA1219588 = D1.1; PRJNA1102327 = the main B3.13 corpus) and
by consensus vs sub-consensus, so the manuscript can separate:
  (a) genotype-defining reference differences (D1.1 consensus calls), from
  (b) genuine within-host sub-consensus events in a 631L/701D (B3.13) background.

Inputs (read-only):
  results/adaptation_site_detections.csv
  data/manifests/v1_20260525_all_hosts.tsv
Outputs:
  results/genotype_stratification_631_701.csv      (per-detection, annotated)
  results/genotype_stratification_summary.json      (counts used in the revised text)

No variant re-calling is performed; this is a re-aggregation of existing per-detection
data. All numbers written to the manuscript revision are taken from the summary here.
"""
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DET = ROOT / "results" / "adaptation_site_detections.csv"
MAN = ROOT / "data" / "manifests" / "v1_20260525_all_hosts.tsv"
OUT_CSV = ROOT / "results" / "genotype_stratification_631_701.csv"
OUT_JSON = ROOT / "results" / "genotype_stratification_summary.json"

# BioProject -> genotype proxy. PRJNA1219588 is the D1.1 genotype project; the main
# USDA corpus PRJNA1102327 is B3.13. Smaller projects are labelled "other".
GENOTYPE = {
    "PRJNA1219588": "D1.1",
    "PRJNA1102327": "B3.13",
}

# Reference (A/cattle/Texas/24-009308-004, B3.13) allele at each site, and whether the
# genotype natively carries the mammalian-adapted allele. This is what determines
# whether a reference-relative difference is a genuine reversion or a genotype artefact.
#   PB2-631: mammalian-adapted = L (M631L). B3.13 carries L; D1.1 is ancestrally M.
#            => 631 L->M in D1.1 is NOT a reversion (D1.1 never had L).
#   PB2-701: mammalian-adapted = N (D701N, forward). B3.13 reference carries D (avian).
#            => 701 D->N in D1.1 consensus is a genotype-level adaptation marker, fixed,
#               not sub-consensus within-host emergence.
SITE_NOTES = {
    "S05": {"gene_pos": "PB2-631", "adapted": "L", "avian": "M", "direction": "reversion"},
    "S02": {"gene_pos": "PB2-701", "adapted": "N", "avian": "D", "direction": "forward"},
}


def classify(row):
    geno = row["genotype"]
    site = row["site_id"]
    subcons = bool(row["is_subconsensus"])
    if site == "S05":  # PB2-631 reversion site
        if geno == "D1.1":
            return "genotype_reference_difference"  # D1.1 ancestrally 631M; not a reversion
        # B3.13 background genuinely carries 631L, so a 631M call here is a candidate event
        return "within_host_reversion_subconsensus" if subcons else "consensus_631M_in_B313_background"
    if site == "S02":  # PB2-701 forward
        if geno == "D1.1":
            return "genotype_consensus_adaptation_marker"  # 701N fixed in D1.1 genotype
        return "within_host_forward_subconsensus" if subcons else "consensus_701N_in_B313_background"
    return "n/a"


def main():
    det = pd.read_csv(DET)
    man = pd.read_csv(MAN, sep="\t")[
        ["run_accession", "bioproject", "host_category", "isolate", "collection_date", "geo_loc_name"]
    ]
    d = det[det.site_id.isin(["S05", "S02"])].merge(
        man, left_on="sample", right_on="run_accession", how="left"
    ).reset_index(drop=True)
    d = d.assign(genotype=d["bioproject"].map(GENOTYPE).fillna("other"))
    d = d.assign(revised_class=d.apply(classify, axis=1))

    keep = [
        "site_id", "gene", "aa_position", "wt_aa", "adapted_aa", "direction",
        "sample", "bioproject", "genotype", "host_category", "isolate",
        "af_mean", "is_subconsensus", "revised_class",
    ]
    d[keep].sort_values(["site_id", "genotype", "af_mean"]).to_csv(OUT_CSV, index=False)

    summary = {"reference": "A/cattle/Texas/24-009308-004 (B3.13); PB2-631 carries adapted L, PB2-701 carries avian D"}
    for site, meta in SITE_NOTES.items():
        s = d[d.site_id == site]
        by_class = s.revised_class.value_counts().to_dict()
        by_geno = s.groupby("genotype").size().to_dict()
        summary[meta["gene_pos"]] = {
            "total_detections": int(len(s)),
            "by_genotype": {k: int(v) for k, v in by_geno.items()},
            "by_revised_class": {k: int(v) for k, v in by_class.items()},
            "subconsensus_in_B313": int(((s.genotype == "B3.13") & (s.is_subconsensus)).sum()),
            "consensus_in_D1.1": int(((s.genotype == "D1.1") & (~s.is_subconsensus)).sum()),
            "notes": meta,
        }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT_CSV.relative_to(ROOT)} and {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
