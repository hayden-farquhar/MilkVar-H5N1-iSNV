"""
Extract per-sample, per-segment coverage statistics from depth files.

For each sample and each of the 8 segments, counts how many positions
have depth >= 100x. Outputs results/segment_coverage_by_sample.parquet.
"""

import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
COVERAGE_DIR = PROJECT_ROOT / "data" / "processed" / "vastai_results" / "coverage"
RESULTS_DIR = PROJECT_ROOT / "results"

DEPTH_THRESHOLD = 100

SEGMENTS = {
    "PP755964.1": ("PB2", 2280),
    "PP755963.1": ("PB1", 2274),
    "PP755962.1": ("PA", 2151),
    "PP755957.1": ("HA", 1704),
    "PP755960.1": ("NP", 1497),
    "PP755959.1": ("NA", 1410),
    "PP755958.1": ("M", 982),
    "PP755961.1": ("NS", 838),
}

SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


def extract_one(depth_file: Path) -> dict:
    sample = depth_file.name.replace(".depth.tsv.gz", "")

    adequate = {gene: 0 for gene in SEGMENT_ORDER}
    total = {gene: 0 for gene in SEGMENT_ORDER}

    with gzip.open(depth_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            depth = int(parts[2])
            gene = SEGMENTS[chrom][0]
            total[gene] += 1
            if depth >= DEPTH_THRESHOLD:
                adequate[gene] += 1

    row = {"sample": sample}
    segments_adequate = 0
    for gene in SEGMENT_ORDER:
        row[f"{gene}_covered"] = adequate[gene]
        row[f"{gene}_total"] = total[gene]
        if adequate[gene] > 0:
            segments_adequate += 1

    row["segments_with_coverage"] = segments_adequate
    row["genome_covered"] = sum(adequate.values())
    row["genome_total"] = sum(total.values())

    return row


def main():
    depth_files = sorted(COVERAGE_DIR.glob("*.depth.tsv.gz"))
    print(f"Found {len(depth_files)} depth files")
    print(f"Computing per-segment coverage (depth >= {DEPTH_THRESHOLD}x)...")

    rows = []
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(extract_one, f): f for f in depth_files}
        done = 0
        for future in as_completed(futures):
            rows.append(future.result())
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(depth_files)}", flush=True)

    df = pd.DataFrame(rows)
    df = df.sort_values("sample").reset_index(drop=True)
    df.to_parquet(RESULTS_DIR / "segment_coverage_by_sample.parquet", index=False)

    print(f"\nWrote {len(df)} samples to results/segment_coverage_by_sample.parquet")

    meets_6seg = (df["segments_with_coverage"] >= 6).sum()
    print(f"Samples with >= 6 segments covered: {meets_6seg}/{len(df)}")

    for gene in SEGMENT_ORDER:
        col = f"{gene}_covered"
        has_any = (df[col] > 0).sum()
        median_frac = (df[col] / df[f"{gene}_total"]).median()
        print(f"  {gene}: {has_any}/{len(df)} with any coverage, "
              f"median {median_frac:.1%} of positions >= {DEPTH_THRESHOLD}x")


if __name__ == "__main__":
    main()
