"""
Concordance filter: intersect iVar and LoFreq variant calls.

A variant is "concordant" when both callers report:
  - Same genomic position
  - Same alternative allele
  - AF >= candidate threshold

Also applies strand-bias filter (Fisher exact test p < threshold → exclude).
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact


def parse_ivar(path: Path) -> pd.DataFrame:
    """Parse iVar variants TSV output."""
    df = pd.read_csv(path, sep="\t")
    # iVar columns: REGION, POS, REF, ALT, REF_DP, REF_RV, ALT_DP, ALT_RV, ...
    # REF_DP = ref forward, REF_RV = ref reverse, ALT_DP = alt forward, ALT_RV = alt reverse
    df = df.rename(columns={
        "REGION": "chrom",
        "POS": "pos",
        "REF": "ref",
        "ALT": "alt",
        "ALT_FREQ": "af",
        "TOTAL_DP": "depth",
        "ALT_DP": "alt_fwd",
        "ALT_RV": "alt_rev",
        "REF_DP": "ref_fwd",
        "REF_RV": "ref_rev",
    })
    # Filter to SNVs only (single nucleotide)
    df = df[df["ref"].str.len() == 1]
    df = df[df["alt"].str.len() == 1]
    df = df[df["alt"] != "+"]  # exclude insertions
    df = df[df["alt"] != "-"]  # exclude deletions
    return df[["chrom", "pos", "ref", "alt", "af", "depth",
               "alt_fwd", "alt_rev", "ref_fwd", "ref_rev"]].copy()


def parse_lofreq(path: Path) -> pd.DataFrame:
    """Parse LoFreq VCF output."""
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, alt, qual, filt, info = fields[:8]
            # Skip indels
            if len(ref) != 1 or len(alt) != 1:
                continue
            # Parse INFO field for AF and DP
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v
            af = float(info_dict.get("AF", 0))
            dp = int(info_dict.get("DP", 0))
            dp4 = info_dict.get("DP4", "0,0,0,0").split(",")
            records.append({
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "af": af,
                "depth": dp,
                "ref_fwd": int(dp4[0]),
                "ref_rev": int(dp4[1]),
                "alt_fwd": int(dp4[2]),
                "alt_rev": int(dp4[3]),
            })
    return pd.DataFrame(records)


def strand_bias_filter(row, p_threshold: float = 0.001) -> bool:
    """Return True if variant PASSES strand-bias filter (not biased)."""
    table = [[row["ref_fwd"], row["ref_rev"]],
             [row["alt_fwd"], row["alt_rev"]]]
    if min(row["alt_fwd"], row["alt_rev"]) == 0 and (row["alt_fwd"] + row["alt_rev"]) < 5:
        return False  # too few alt reads to assess
    try:
        _, p = fisher_exact(table)
        return p >= p_threshold
    except ValueError:
        return True


def main(ivar_path, lofreq_path, output_path, min_freq=0.01, strand_bias_p=0.001):
    ivar_df = parse_ivar(Path(ivar_path))
    lofreq_df = parse_lofreq(Path(lofreq_path))

    if ivar_df.empty or lofreq_df.empty:
        pd.DataFrame(columns=[
            "chrom", "pos", "ref", "alt", "af_ivar", "af_lofreq",
            "depth_ivar", "depth_lofreq", "concordant", "passes_strand_bias"
        ]).to_csv(output_path, sep="\t", index=False)
        return

    # Merge on position + alt allele (concordance definition)
    merged = ivar_df.merge(
        lofreq_df,
        on=["chrom", "pos", "ref", "alt"],
        suffixes=("_ivar", "_lofreq"),
        how="inner",
    )

    # Both must be at >= min_freq
    merged = merged[
        (merged["af_ivar"] >= min_freq) & (merged["af_lofreq"] >= min_freq)
    ]

    # Apply strand-bias filter using iVar strand counts
    merged["passes_strand_bias"] = merged.apply(
        lambda row: strand_bias_filter(row, strand_bias_p), axis=1
    )

    # Mean AF across callers
    merged["af_mean"] = (merged["af_ivar"] + merged["af_lofreq"]) / 2

    # Output columns
    out = merged[[
        "chrom", "pos", "ref", "alt",
        "af_ivar", "af_lofreq", "af_mean",
        "depth_ivar", "depth_lofreq",
        "alt_fwd_ivar", "alt_rev_ivar",
        "passes_strand_bias",
    ]].copy()
    out["concordant"] = True

    out.to_csv(output_path, sep="\t", index=False)
    print(f"  {Path(output_path).stem}: {len(out)} concordant SNVs "
          f"({out['passes_strand_bias'].sum()} pass strand-bias filter)")


if __name__ == "__main__":
    import snakemake
    main(
        snakemake.input.ivar,
        snakemake.input.lofreq,
        snakemake.output.concordant,
        min_freq=snakemake.params.min_freq,
        strand_bias_p=snakemake.params.strand_bias_p,
    )
