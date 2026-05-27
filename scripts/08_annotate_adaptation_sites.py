"""
Annotate concordant iSNVs against the Tier 1 mammalian-adaptation site panel.

Maps genomic coordinates to canonical amino acid positions and flags
variants at pre-specified adaptation sites.
"""

import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


TIER1_SITES = [
    {"id": "S01", "gene": "PB2", "segment": 1, "aa_pos": 627, "wt": "E", "adapted": "K"},
    {"id": "S02", "gene": "PB2", "segment": 1, "aa_pos": 701, "wt": "D", "adapted": "N"},
    {"id": "S03", "gene": "PB2", "segment": 1, "aa_pos": 591, "wt": "Q", "adapted": "K"},
    {"id": "S04", "gene": "PB2", "segment": 1, "aa_pos": 271, "wt": "T", "adapted": "A"},
    {"id": "S05", "gene": "PB2", "segment": 1, "aa_pos": 631, "wt": "M", "adapted": "L"},
    {"id": "S06", "gene": "PA", "segment": 3, "aa_pos": 497, "wt": "K", "adapted": "R"},
    {"id": "S07", "gene": "HA", "segment": 4, "aa_pos": None, "wt": "Q", "adapted": "L",
     "h3_pos": 226, "note": "H5 position determined by alignment"},
    {"id": "S08", "gene": "HA", "segment": 4, "aa_pos": None, "wt": "G", "adapted": "S",
     "h3_pos": 228, "note": "H5 position determined by alignment"},
    {"id": "S09", "gene": "PB1-F2", "segment": 2, "aa_pos": 66, "wt": "N", "adapted": "S"},
    {"id": "S10", "gene": "NS1", "segment": 8, "aa_pos": 92, "wt": "D", "adapted": "E"},
    {"id": "S11", "gene": "M2", "segment": 7, "aa_pos": 31, "wt": "S", "adapted": "N"},
]

# Segment accession to segment number mapping
SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


def load_reference(fasta_path: Path) -> dict:
    """Load reference genome and map accessions to segments."""
    seqs = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seqs[rec.id] = str(rec.seq)
    return seqs


def codon_at_position(seq: str, aa_pos: int, frame_start: int = 0) -> tuple:
    """Get the codon and nucleotide positions for a given amino acid position.

    Returns (codon_str, nt_start_0based, nt_end_0based).
    """
    nt_start = frame_start + (aa_pos - 1) * 3
    nt_end = nt_start + 3
    if nt_end > len(seq):
        return None, None, None
    codon = seq[nt_start:nt_end]
    return codon, nt_start, nt_end


def find_orf_start(seq: str, gene: str) -> int:
    """Find the start of the main ORF for a given gene/segment.

    For influenza segments, the main ORF starts at the first ATG.
    This is a simplified approach; a proper annotation would use
    GenBank feature tables.
    """
    # Most influenza ORFs start near position 1-30
    # Find first ATG in the first 100 bases
    for i in range(0, min(100, len(seq) - 2)):
        if seq[i:i+3] == "ATG":
            return i
    return 0


def build_site_coordinate_map(ref_seqs: dict, ref_metadata_path: Path) -> list:
    """Build the coordinate map from adaptation sites to genomic positions."""
    meta = json.loads(ref_metadata_path.read_text())
    segment_accs = meta["segments"]

    # Map segment names to accessions
    acc_to_segment = {v: k for k, v in segment_accs.items()}

    coordinate_map = []

    for site in TIER1_SITES:
        gene = site["gene"]
        aa_pos = site["aa_pos"]

        if aa_pos is None:
            # HA sites with H3 numbering — need alignment-based mapping
            # For now, mark as requiring manual verification
            coordinate_map.append({
                **site,
                "nt_positions": None,
                "status": "requires_alignment_verification"
            })
            continue

        # Find the right segment
        if gene == "PB2":
            seg_name = "PB2"
        elif gene == "PB1" or gene == "PB1-F2":
            seg_name = "PB1"
        elif gene == "PA":
            seg_name = "PA"
        elif gene == "HA":
            seg_name = "HA"
        elif gene == "NP":
            seg_name = "NP"
        elif gene == "NA":
            seg_name = "NA"
        elif gene in ("M1", "M2"):
            seg_name = "M"
        elif gene in ("NS1", "NS2", "NEP"):
            seg_name = "NS"
        else:
            seg_name = gene

        acc = segment_accs.get(seg_name)
        if acc is None or acc not in ref_seqs:
            coordinate_map.append({**site, "nt_positions": None, "status": "segment_not_found"})
            continue

        seq = ref_seqs[acc]
        orf_start = find_orf_start(seq, gene)

        # Special handling for M2 (overlapping reading frame on segment 7)
        # M2 uses a spliced mRNA; the M2 ORF starts at a different position
        if gene == "M2":
            # M2 shares the first 26 aa with M1, then continues from the splice junction
            # For position 31, we're in the shared N-terminal region
            # M2 position 31 corresponds to M1 position 31 in the shared reading frame
            pass  # Use same ORF start as M1

        # Special handling for PB1-F2 (+1 frameshift on segment 2)
        if gene == "PB1-F2":
            # PB1-F2 is in the +1 reading frame relative to PB1
            # Its start codon is ~95-120 nt downstream of PB1 start in +1 frame
            # Find the first ATG in +1 frame after position ~90
            for i in range(orf_start + 91, min(orf_start + 200, len(seq) - 2)):
                if seq[i:i+3] == "ATG":
                    orf_start = i
                    break

        codon, nt_start, nt_end = codon_at_position(seq, aa_pos, orf_start)

        if codon is None:
            coordinate_map.append({**site, "nt_positions": None, "status": "out_of_range"})
            continue

        # Translate to verify wildtype amino acid
        try:
            translated_aa = str(Seq(codon).translate())
        except Exception:
            translated_aa = "?"

        wt_match = translated_aa == site["wt"]

        coordinate_map.append({
            **site,
            "segment_acc": acc,
            "orf_start": orf_start,
            "nt_positions": list(range(nt_start, nt_end)),
            "codon": codon,
            "translated_aa": translated_aa,
            "wt_verified": wt_match,
            "status": "verified" if wt_match else "wt_mismatch"
        })

    return coordinate_map


def annotate_variants(variants_path: Path, coordinate_map: list) -> pd.DataFrame:
    """Annotate concordant variants against the adaptation-site map."""
    if not variants_path.exists():
        return pd.DataFrame()

    variants = pd.read_csv(variants_path, sep="\t")
    if variants.empty:
        return pd.DataFrame()

    annotations = []
    for site in coordinate_map:
        if site.get("nt_positions") is None:
            continue

        acc = site.get("segment_acc", "")
        for nt_pos in site["nt_positions"]:
            # Find variants at this position (0-based in depth file, 1-based in VCF/iVar)
            matches = variants[
                (variants["chrom"] == acc) &
                (variants["pos"] == nt_pos + 1)  # iVar uses 1-based positions
            ]
            for _, var in matches.iterrows():
                annotations.append({
                    "site_id": site["id"],
                    "gene": site["gene"],
                    "aa_position": site["aa_pos"],
                    "wt_aa": site["wt"],
                    "adapted_aa": site["adapted"],
                    "chrom": var["chrom"],
                    "nt_position": var["pos"],
                    "ref_nt": var["ref"],
                    "alt_nt": var["alt"],
                    "af_mean": var.get("af_mean", var.get("af_ivar", 0)),
                    "depth": var.get("depth_ivar", 0),
                    "passes_strand_bias": var.get("passes_strand_bias", True),
                })

    return pd.DataFrame(annotations)


if __name__ == "__main__":
    project_root = Path(__file__).resolve().parent.parent
    ref_path = project_root / "data" / "reference" / "h5n1_b3.13_cattle_texas_reference.fasta"
    meta_path = project_root / "data" / "reference" / "reference_metadata.json"

    ref_seqs = load_reference(ref_path)
    coord_map = build_site_coordinate_map(ref_seqs, meta_path)

    print("Adaptation-site coordinate map:")
    for site in coord_map:
        status = site.get("status", "unknown")
        pos_str = site.get("nt_positions", "N/A")
        if pos_str and isinstance(pos_str, list):
            pos_str = f"nt {pos_str[0]+1}-{pos_str[-1]+1}"
        print(f"  {site['id']} ({site['gene']} {site.get('aa_pos','?')}): "
              f"{pos_str} | codon={site.get('codon','?')} | "
              f"translated={site.get('translated_aa','?')} | {status}")
