#!/usr/bin/env python3
"""Merge Vast.ai results into corpus_variants.parquet and corpus_coverage.parquet.

Run after rsyncing results/ back from Vast.ai:
    python3 scripts/merge_results.py --results data/processed/vastai_results --manifest data/manifests/v1_20260525.tsv
"""
import argparse
import gzip
import json
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from tqdm import tqdm

MIN_FREQ = 0.01


def parse_ivar(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep='\t')
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return pd.DataFrame()
    snvs = df[
        (df['ALT_FREQ'] >= MIN_FREQ)
        & (df['REF'].str.len() == 1)
        & (df['ALT'].str.len() == 1)
        & (~df['ALT'].isin(['+', '-']))
    ].copy()
    snvs['key'] = snvs['REGION'] + ':' + snvs['POS'].astype(str) + ':' + snvs['ALT']
    return snvs


def parse_lofreq(path: Path) -> pd.DataFrame:
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8 or len(parts[3]) != 1 or len(parts[4]) != 1:
                continue
            info = dict(x.split('=', 1) for x in parts[7].split(';') if '=' in x)
            af = float(info.get('AF', 0))
            dp = int(info.get('DP', 0))
            if af >= MIN_FREQ:
                records.append({
                    'chrom': parts[0], 'pos': int(parts[1]), 'alt': parts[4],
                    'af_lofreq': af, 'depth_lofreq': dp,
                    'key': f'{parts[0]}:{parts[1]}:{parts[4]}'
                })
    return pd.DataFrame(records) if records else pd.DataFrame()


def concordance(accession: str, variants_dir: Path) -> pd.DataFrame:
    ivar_path = variants_dir / f'{accession}.ivar.tsv'
    lofreq_path = variants_dir / f'{accession}.lofreq.vcf'
    if not ivar_path.exists() or not lofreq_path.exists():
        return pd.DataFrame()

    ivar_snvs = parse_ivar(ivar_path)
    lofreq_df = parse_lofreq(lofreq_path)
    if ivar_snvs.empty or lofreq_df.empty:
        return pd.DataFrame()

    merged = ivar_snvs.merge(lofreq_df, on='key', how='inner')
    if merged.empty:
        return pd.DataFrame()

    def strand_bias_ok(row):
        try:
            _, p = fisher_exact([
                [row.get('REF_DP', 0), row.get('REF_RV', 0)],
                [row.get('ALT_DP', 0), row.get('ALT_RV', 0)]
            ])
            return p >= 0.001
        except Exception:
            return True

    merged['passes_strand_bias'] = merged.apply(strand_bias_ok, axis=1)
    merged['af_mean'] = (merged['ALT_FREQ'] + merged['af_lofreq']) / 2
    merged['sample'] = accession

    return merged[['sample', 'REGION', 'POS', 'REF', 'ALT',
                    'ALT_FREQ', 'af_lofreq', 'af_mean',
                    'TOTAL_DP', 'depth_lofreq',
                    'passes_strand_bias']].rename(columns={
        'REGION': 'chrom', 'POS': 'pos', 'REF': 'ref', 'ALT': 'alt',
        'ALT_FREQ': 'af_ivar', 'TOTAL_DP': 'depth_ivar'
    })


def median_depth(depth_path: Path) -> int:
    depths = []
    opener = gzip.open if str(depth_path).endswith('.gz') else open
    with opener(depth_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))
    if not depths:
        return 0
    depths.sort()
    return depths[len(depths) // 2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', required=True, help='Path to rsynced results/')
    parser.add_argument('--manifest', required=True, help='Path to manifest TSV')
    parser.add_argument('--output', default='results', help='Output directory')
    args = parser.parse_args()

    results = Path(args.results)
    variants_dir = results / 'variants'
    coverage_dir = results / 'coverage'
    qc_dir = results / 'qc'
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep='\t')
    completed_file = results / 'completed.txt'
    if completed_file.exists():
        completed = set(completed_file.read_text().strip().split('\n'))
    else:
        completed = {p.stem.replace('.ivar', '') for p in variants_dir.glob('*.ivar.tsv')}

    print(f'Completed samples: {len(completed)}')
    print(f'Manifest samples:  {len(manifest)}')

    # Coverage stats
    print('\nComputing coverage stats...')
    stats_rows = []
    for acc in tqdm(completed, desc='Coverage'):
        depth_path = coverage_dir / f'{acc}.depth.tsv.gz'
        ivar_path = variants_dir / f'{acc}.ivar.tsv'
        row_data = manifest[manifest['run_accession'] == acc]

        md = median_depth(depth_path) if depth_path.exists() else 0
        n_ivar = 0
        if ivar_path.exists():
            try:
                idf = pd.read_csv(ivar_path, sep='\t')
                n_ivar = len(idf[idf['ALT_FREQ'] >= MIN_FREQ]) if not idf.empty else 0
            except Exception:
                pass

        stats_rows.append({
            'run_accession': acc,
            'host_category': row_data['host_category'].iloc[0] if len(row_data) > 0 else 'unknown',
            'platform': row_data['platform'].iloc[0] if len(row_data) > 0 else 'unknown',
            'median_depth': md,
            'ivar_snvs': n_ivar,
        })

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_parquet(output / 'corpus_coverage.parquet', index=False)
    print(f'\nDepth summary:')
    print(stats_df.groupby('host_category')['median_depth'].describe().to_string())
    print(f'\nBelow 100x: {(stats_df["median_depth"] < 100).sum()}')

    # Concordance
    print('\nComputing concordance...')
    all_concordant = []
    for acc in tqdm(completed, desc='Concordance'):
        conc = concordance(acc, variants_dir)
        if not conc.empty:
            all_concordant.append(conc)

    if all_concordant:
        corpus = pd.concat(all_concordant, ignore_index=True)
        corpus.to_parquet(output / 'corpus_variants.parquet', index=False)
        print(f'\nCorpus concordant variants: {len(corpus):,}')
        print(f'Samples with concordant calls: {corpus["sample"].nunique()}')
        print(f'Passing strand bias: {corpus["passes_strand_bias"].sum():,}')
    else:
        print('No concordant variants found.')

    # Failed samples
    failed_file = results / 'failed.txt'
    if failed_file.exists():
        failed = failed_file.read_text().strip().split('\n')
        failed = [f for f in failed if f]
        print(f'\nFailed samples: {len(failed)}')
        if failed:
            (output / 'failed_samples.txt').write_text('\n'.join(failed))


if __name__ == '__main__':
    main()
