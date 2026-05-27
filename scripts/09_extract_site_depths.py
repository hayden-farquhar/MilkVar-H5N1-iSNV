"""
Extract per-sample depth at all 11 Tier 1 adaptation-site codon positions
from the per-sample samtools depth files.

Reads 4,559 gzipped depth TSVs, extracts 33 positions (11 sites × 3 nt),
computes min codon depth per site per sample, and writes:
  results/site_depth_by_sample.parquet
"""

import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
COVERAGE_DIR = PROJECT_ROOT / "data" / "processed" / "vastai_results" / "coverage"
RESULTS_DIR = PROJECT_ROOT / "results"

SITE_POSITIONS = {
    "S01": ("PP755964.1", [1879, 1880, 1881]),
    "S02": ("PP755964.1", [2101, 2102, 2103]),
    "S03": ("PP755964.1", [1771, 1772, 1773]),
    "S04": ("PP755964.1", [811, 812, 813]),
    "S05": ("PP755964.1", [1891, 1892, 1893]),
    "S06": ("PP755962.1", [1489, 1490, 1491]),
    "S07": ("PP755957.1", [712, 713, 714]),
    "S08": ("PP755957.1", [718, 719, 720]),
    "S09": ("PP755963.1", [290, 291, 292]),
    "S10": ("PP755961.1", [274, 275, 276]),
    "S11": ("PP755958.1", [779, 780, 781]),
}

TARGET_KEYS = set()
POS_TO_SITES = {}
for site_id, (chrom, positions) in SITE_POSITIONS.items():
    for pos in positions:
        key = (chrom, pos)
        TARGET_KEYS.add(key)
        POS_TO_SITES.setdefault(key, []).append(site_id)


def extract_one(depth_file: Path) -> dict:
    sample = depth_file.name.replace(".depth.tsv.gz", "")
    site_depths = {sid: [] for sid in SITE_POSITIONS}

    with gzip.open(depth_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            key = (chrom, pos)
            if key in TARGET_KEYS:
                depth = int(parts[2])
                for sid in POS_TO_SITES[key]:
                    site_depths[sid].append(depth)

    row = {"sample": sample}
    for sid in SITE_POSITIONS:
        depths = site_depths[sid]
        row[f"{sid}_min_depth"] = min(depths) if len(depths) == 3 else 0
    return row


def main():
    depth_files = sorted(COVERAGE_DIR.glob("*.depth.tsv.gz"))
    print(f"Found {len(depth_files)} depth files")
    print(f"Extracting depth at {len(TARGET_KEYS)} positions across 11 sites...")

    rows = []
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(extract_one, f): f for f in depth_files}
        done = 0
        for future in as_completed(futures):
            rows.append(future.result())
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(depth_files)}")

    df = pd.DataFrame(rows)
    df = df.sort_values("sample").reset_index(drop=True)
    df.to_parquet(RESULTS_DIR / "site_depth_by_sample.parquet", index=False)

    print(f"\nWrote {len(df)} samples to results/site_depth_by_sample.parquet")
    for sid in SITE_POSITIONS:
        col = f"{sid}_min_depth"
        adequate = (df[col] >= 100).sum()
        print(f"  {sid}: {adequate}/{len(df)} samples with >=100x depth")


if __name__ == "__main__":
    main()
