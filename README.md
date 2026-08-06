# No receptor-binding domain adaptation detected in within-host H5N1 surveillance of 4,559 US dairy outbreak sequences

Code repository for: **No receptor-binding domain adaptation detected in within-host H5N1 surveillance of 4,559 US dairy outbreak sequences**

Hayden Farquhar MBBS MPHTM, Independent researcher, Finley, NSW, Australia

Pre-registration: [10.17605/OSF.IO/WY9FA](https://doi.org/10.17605/OSF.IO/WY9FA)

## Overview

This repository contains the analysis code for a pre-registered, corpus-wide intrahost single-nucleotide variant (iSNV) analysis of all publicly available H5N1 cattle, feline-spillover, and retail-milk sequences on the NCBI SRA (4,559 samples). The analysis surveys an 11-site mammalian-adaptation panel, compares within-host diversity across host categories, and provides a sub-consensus surveillance baseline for the US dairy H5N1 outbreak.

## Data Sources

| Source | URL | Access | Notes |
|--------|-----|--------|-------|
| NCBI SRA | https://www.ncbi.nlm.nih.gov/sra | Public | Raw FASTQ reads; 7 BioProjects listed in `data/raw/manifest.tsv` |
| NCBI GenBank | https://www.ncbi.nlm.nih.gov/genbank | Public | Reference genome PP755957-PP755964 |

Raw sequencing reads (>500 GB total) are not redistributed. They can be downloaded from the NCBI SRA using the accessions in `data/raw/manifest.tsv` via `prefetch` and `fasterq-dump` (sra-tools) or ENA download.

The processed variant tables (`data/processed/`) are included in this repository (13 MB total) to allow reproduction of all analyses from Phase 5 onwards without re-running the full variant-calling pipeline.

## Requirements

**For analysis scripts (Phase 5+):** Python 3.11+

```bash
pip install -r requirements.txt
```

**For the full pipeline (Phases 1-4):** additionally requires BWA, minimap2, samtools, iVar, LoFreq, fastp, and sra-tools. See `data/raw/README.md` for installation instructions.

## Reproduction

### Quick path (from processed data)

The processed variant tables are included. To reproduce all reported results and figures:

```bash
# Step 1: Adaptation-site depth extraction
python scripts/09_extract_site_depths.py

# Step 2: Segment-level coverage extraction
python scripts/10_extract_segment_coverage.py

# Step 3: Adaptation-site prevalence analysis (Table 2, Figure 2)
python scripts/11_adaptation_analysis.py

# Step 4: Diversity comparisons (Figure 3)
python scripts/12_diversity_analysis.py

# Step 5: Generate figures
python scripts/13_figure_heatmap.py
python scripts/14_figure_corpus_qc.py
python scripts/15_figure_threshold.py
python scripts/16_figure_variant_density.py
```

Estimated runtime: ~15 minutes (Steps 1-2 are I/O bound; Steps 3-5 are fast).

Note: Steps 1-2 require the per-sample depth files (`data/processed/vastai_results/coverage/*.depth.tsv.gz`), which are 46 MB compressed and are NOT included in this repository due to size. The pre-extracted outputs (`data/processed/site_depth_by_sample.parquet` and `data/processed/segment_coverage_by_sample.parquet`) are included, so Steps 3-5 can run without Steps 1-2.

### Full pipeline (from raw SRA reads)

To reproduce the entire analysis from raw data:

```bash
# Phase 1: Build SRA manifest
python scripts/01_build_manifest.py

# Phase 3: Process each sample (requires BWA/minimap2/samtools/iVar/LoFreq)
# This processes 4,559 samples. On a 72-core machine it takes ~25 hours.
bash scripts/03_run_corpus.sh

# Phase 3: Merge per-sample results
python scripts/05_merge_results.py --results data/processed/vastai_results --manifest data/raw/manifest.tsv

# Phase 4: Spike-in validation (requires processed BAMs)
# See scripts/06_spikein_validation.py for setup

# Phase 5+: Continue with Quick path above
```

## Script Descriptions

| Script | Phase | Description | Inputs | Outputs |
|--------|-------|-------------|--------|---------|
| `01_build_manifest.py` | 1 | Query NCBI Entrez for H5N1 cattle/feline/milk SRA runs | NCBI API | `data/raw/manifest.tsv` |
| `02_process_sample.sh` | 3 | Per-sample pipeline: fetch, QC, align, depth, iVar, LoFreq | SRA accession | Per-sample variant calls + depth |
| `03_run_corpus.sh` | 3 | Orchestrate corpus-scale processing via GNU parallel | Manifest | All per-sample results |
| `04_concordance_filter.py` | 3 | Dual-caller concordance intersection | iVar + LoFreq calls | Concordant variants |
| `05_merge_results.py` | 3 | Merge per-sample results into corpus parquets | Per-sample results | `corpus_variants.parquet`, `corpus_coverage.parquet` |
| `06_spikein_validation.py` | 4 | V3 synthetic spike-in detection-limit validation | Two high-coverage BAMs | Spike-in sensitivity table |
| `07_consensus_validation.sh` | 4 | V1 consensus recovery validation | Selected SRA runs | Consensus comparison |
| `08_annotate_adaptation_sites.py` | 5 | Map genomic coordinates to Tier 1 adaptation sites | Reference genome | Site coordinate definitions |
| `09_extract_site_depths.py` | 5 | Extract per-sample depth at 11 adaptation-site codons | 4,559 depth files | `site_depth_by_sample.parquet` |
| `10_extract_segment_coverage.py` | 6 | Extract per-segment coverage statistics | 4,559 depth files | `segment_coverage_by_sample.parquet` |
| `11_adaptation_analysis.py` | 5 | Adaptation-site prevalence with CIs and Bonferroni correction | Variants + depths | Tables 2-3, prevalence CSVs |
| `12_diversity_analysis.py` | 6 | Nucleotide diversity, entropy, iSNV counts; Kruskal-Wallis tests | Variants + coverage | Figure 3, diversity parquet, test results |
| `13_figure_heatmap.py` | 5 | Figure 2: adaptation-site allele frequency heatmap | Detections + depths | `figure1_adaptation_heatmap.{png,tiff}` |
| `14_figure_corpus_qc.py` | — | Figure 4: corpus composition and QC summary | Coverage + manifest | `figure3_corpus_qc.{png,tiff}` |
| `15_figure_threshold.py` | 4 | Figure 1: threshold determination criteria plot | Phase 4 results; C2 computed from variant table | `figure4_threshold.{png,tiff}` |
| `16_figure_variant_density.py` | 6 | Suppl. Figure S1: genome-wide variant density profile | Sub-consensus variants | `figure_s1_variant_density.{png,tiff}` |
| `17_genotype_stratification.py` | 5 | Stratify PB2-631/701 detections by genotype (B3.13 vs D1.1) | Detections + BioProject | `genotype_stratification_*.{csv,json}` |
| `18_rebuild_corpus_variants.py` | 3 | Rebuild the concordant variant table under the corrected strand-bias filter | Per-caller outputs | `corpus_variants*.parquet`, `strandfix_delta.json` |
| `19_sa2_sa5_sensitivity.py` | — | SA2-SA5: single-caller, stringent depth, platform, library strategy | Per-caller outputs + variants | `sa2_*.csv`, `sa2_sa5_summary.json` |
| `20_sensitivity_analyses.py` | — | SA1, SA6, non-synonymous sensitivity, mixed-infection screen | Variants + diversity | `sa1_*.csv`, `sa6_*.json`, `sensitivity_all_nonsyn_*.csv`, `mixed_infection_flagged.csv` |
| `21_tier2_exploratory_screen.py` | 5 | Pre-registered Tier 2 exploratory site screen | Variants + reference | `tier2_screen_*.{csv,json}` |
| `22_build_table_s5.py` | — | Render Supplementary Table S5 | SA2-SA5 outputs | `table_s5_sa2_sa5.md` |
| `23_build_table_s6_qc.py` | — | Render Supplementary Table S6; parameters parsed from the pipeline definitions | `config.yaml`, `rules/*.smk`, coverage | `table_s6_qc_metrics.{csv,md}` |
| `_pipeline.py` | — | Loader giving scripts 17+ importable access to the numbered modules | — | — |

Scripts 17 onward reuse the site definitions, codon classification, and caller
parsers from the earlier stages via `_pipeline.py`, rather than restating them.
The numeric prefixes document run order but are not valid Python module names,
which is what the loader works around.

### Note on version 2.0.0

Release 2.0.0 corrects a defect in the strand-bias filter present in 1.0.0.
iVar reports `REF_DP`/`ALT_DP` as total depth per allele and `REF_RV`/`ALT_RV`
as the reverse-strand portion; both `04_concordance_filter.py` and
`05_merge_results.py` passed the totals into Fisher's exact test as though they
were forward-strand counts, which left the filter too permissive. Correcting it
reclassified 16,067 of 289,054 concordant calls (5.6%), all in the same
direction. `11_adaptation_analysis.py` additionally had the zero-count branch of
its confidence-interval calculation using the Bonferroni alpha rather than the
95% level used elsewhere. Results derived from 1.0.0 should not be used; see
`strandfix_delta.json` for the full delta.

## Outputs

| File | Manuscript reference |
|------|---------------------|
| `outputs/tables/adaptation_site_prevalence.csv` | Table 2 |
| `outputs/tables/adaptation_cooccurrence.csv` | Table 3 (co-occurrence data) |
| `outputs/tables/diversity_tests.json` | Results section 4.6 |
| `outputs/tables/sa1_threshold_sensitivity.csv` | Supplementary Table S2 |
| `outputs/tables/sensitivity_all_nonsyn_summary.csv` | Supplementary Table S3 |
| `outputs/tables/sensitivity_all_nonsyn_detections.csv` | Supplementary Table S4 |
| `outputs/figures/figure1_adaptation_heatmap.png` | Figure 2 (preprint numbering) |
| `outputs/figures/figure2_diversity.png` | Figure 3 (preprint numbering) |
| `outputs/figures/figure3_corpus_qc.png` | Figure 4 (preprint numbering) |
| `outputs/figures/figure4_threshold.png` | Figure 1 (preprint numbering) |
| `outputs/figures/figure_s1_variant_density.png` | Supplementary Figure S1 |

## Citation

If you use this code, please cite:

```
Farquhar H. No receptor-binding domain adaptation detected in within-host H5N1
surveillance of 4,559 US dairy outbreak sequences. Zenodo preprint, 2026.
```

## License

Code: MIT License. Data and documentation: CC-BY 4.0.
