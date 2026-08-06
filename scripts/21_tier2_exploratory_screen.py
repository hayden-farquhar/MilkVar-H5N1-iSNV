"""
Tier 2 exploratory site screen (pre-registration Section 9.1).

The pre-registration catalogued a set of Tier 2 sites of biological interest that
carried no hypothesis test and were to be reported descriptively if detected at
AF >= threshold. This screen was not executed in the original analysis; this
script performs it against the corrected concordant variant table.

Coordinate mapping notes:
  * The reference segments are CDS-aligned from nt 1, verified against the Tier 1
    panel (HA full-length 238 -> Q, 240 -> G as expected). A residue N therefore
    occupies nt 3N-2 .. 3N.
  * PA-X shares its first ~191 residues with PA (the +1 ribosomal frameshift occurs
    near codon 191), so PA-X residue 42 is identical to PA residue 42.
  * The pre-registration specifies HA residues 158/160/182/192 in "H5 numbering"
    without stating whether this is full-length or mature (signal-peptide-trimmed)
    numbering, and the reference translation does not disambiguate it. Both
    interpretations are therefore screened and reported separately.
  * NA stalk truncations (residues 49-68) are deletions. The analysis pipeline
    retains single-nucleotide variants only, so these are out of scope and are
    reported as not assessable rather than as absent.

Any non-synonymous variant at a Tier 2 codon is reported (not a single
pre-specified substitution), consistent with the exploratory framing.

Outputs:
  results/tier2_screen_detections.csv
  results/tier2_screen_summary.csv
  results/tier2_screen_summary.json
"""

import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
REFERENCE = PROJECT_ROOT / "data" / "reference" / "h5n1_b3.13_cattle_texas_reference.fasta"

AF_THRESHOLD = 0.03
DEPTH_THRESHOLD = 100
H5_SIGNAL_PEPTIDE = 16

SEGMENTS = {"HA": "PP755957.1", "NP": "PP755960.1", "PA": "PP755962.1"}

# (label, gene, accession, residue position in the CDS, provenance note)
TIER2_TARGETS = [
    ("HA-158 (full-length numbering)", "HA", "HA", 158, "H5 numbering read as full-length"),
    ("HA-160 (full-length numbering)", "HA", "HA", 160, "H5 numbering read as full-length"),
    ("HA-182 (full-length numbering)", "HA", "HA", 182, "H5 numbering read as full-length"),
    ("HA-192 (full-length numbering)", "HA", "HA", 192, "H5 numbering read as full-length"),
    ("HA-158 (mature numbering)", "HA", "HA", 158 + H5_SIGNAL_PEPTIDE, "H5 numbering read as mature protein"),
    ("HA-160 (mature numbering)", "HA", "HA", 160 + H5_SIGNAL_PEPTIDE, "H5 numbering read as mature protein"),
    ("HA-182 (mature numbering)", "HA", "HA", 182 + H5_SIGNAL_PEPTIDE, "H5 numbering read as mature protein"),
    ("HA-192 (mature numbering)", "HA", "HA", 192 + H5_SIGNAL_PEPTIDE, "H5 numbering read as mature protein"),
    ("PA-X-42", "PA-X", "PA", 42, "PA-X shares residues 1-191 with PA"),
    ("NP-319", "NP", "NP", 319, "standard ORF"),
    ("NP-357", "NP", "NP", 357, "standard ORF"),
]

NOT_ASSESSABLE = [{
    "site": "NA stalk residues 49-68",
    "gene": "NA",
    "reason": "Stalk truncations are deletions; the pipeline retains SNVs only "
              "(indels excluded at the concordance filter). Not assessable — "
              "absence of a result here is not evidence of absence.",
}]


def load_reference():
    return {r.id: str(r.seq) for r in SeqIO.parse(REFERENCE, "fasta")}


def screen():
    seqs = load_reference()
    variants = pd.read_parquet(RESULTS_DIR / "corpus_variants.parquet")
    coverage = pd.read_parquet(RESULTS_DIR / "corpus_coverage.parquet").rename(
        columns={"run_accession": "sample"}
    )
    host = dict(zip(coverage["sample"], coverage["host_category"]))

    detections, summary = [], []

    for label, gene, seg_key, aa_pos, note in TIER2_TARGETS:
        acc = SEGMENTS[seg_key]
        seq = seqs[acc]
        nt_positions = [3 * aa_pos - 2, 3 * aa_pos - 1, 3 * aa_pos]
        ref_codon = seq[nt_positions[0] - 1: nt_positions[2]]
        ref_aa = str(Seq(ref_codon).translate())

        hits = variants[
            (variants["chrom"] == acc)
            & (variants["pos"].isin(nt_positions))
            & (variants["af_mean"] >= AF_THRESHOLD)
            & (variants["passes_strand_bias"])
        ]

        n_nonsyn = 0
        for v in hits.itertuples(index=False):
            idx = nt_positions.index(v.pos)
            mutant = list(ref_codon)
            mutant[idx] = v.alt
            mut_aa = str(Seq("".join(mutant)).translate())
            if mut_aa == ref_aa:
                continue
            n_nonsyn += 1
            detections.append({
                "site": label, "gene": gene, "aa_position": aa_pos,
                "ref_aa": ref_aa, "alt_aa": mut_aa,
                "sample": v.sample, "host": host.get(v.sample, "unknown"),
                "nt_position": v.pos, "ref_nt": v.ref, "alt_nt": v.alt,
                "af_mean": v.af_mean,
                "sub_consensus": bool(v.af_mean < 0.50),
                "numbering_note": note,
            })

        summary.append({
            "site": label, "gene": gene, "aa_position_in_cds": aa_pos,
            "reference_codon": ref_codon, "reference_aa": ref_aa,
            "n_nonsynonymous_detections": n_nonsyn,
            "numbering_note": note,
        })

    return pd.DataFrame(detections), pd.DataFrame(summary)


def main():
    det, summ = screen()
    det.to_csv(RESULTS_DIR / "tier2_screen_detections.csv", index=False)
    summ.to_csv(RESULTS_DIR / "tier2_screen_summary.csv", index=False)

    payload = {
        "af_threshold": AF_THRESHOLD,
        "depth_threshold": DEPTH_THRESHOLD,
        "total_nonsynonymous_detections": int(len(det)),
        "sites_screened": int(len(summ)),
        "not_assessable": NOT_ASSESSABLE,
    }
    (RESULTS_DIR / "tier2_screen_summary.json").write_text(json.dumps(payload, indent=2))

    print(summ.to_string(index=False))
    print()
    if det.empty:
        print("No non-synonymous Tier 2 detections at AF >= 3%.")
    else:
        print(det.to_string(index=False))
    print(f"\nNot assessable: {NOT_ASSESSABLE[0]['site']} — {NOT_ASSESSABLE[0]['reason']}")


if __name__ == "__main__":
    main()
