# Data Dictionary

## Raw data

### `data/raw/manifest.tsv`

Corpus manifest (4,571 SRA runs). Tab-separated.

| Variable | Type | Description |
|----------|------|-------------|
| `run_accession` | string | NCBI SRA run accession (e.g., SRR28752447) |
| `total_bases` | int | Total sequenced bases |
| `total_spots` | int | Total read pairs/spots |
| `isolate` | string | Isolate name from SRA metadata |
| `collection_date` | string | Collection date (YYYY-MM-DD or YYYY) |
| `geo_loc_name` | string | Geographic location (e.g., "USA: Texas") |
| `host` | string | Host organism (e.g., "Bos taurus", "Felis catus") |
| `biosample` | string | NCBI BioSample accession |
| `library_strategy` | string | WGS or AMPLICON |
| `library_selection` | string | Library selection method |
| `platform` | string | ILLUMINA or OXFORD_NANOPORE |
| `instrument` | string | Sequencing instrument model |
| `bioproject` | string | NCBI BioProject accession |
| `host_category` | string | Assigned host category: cattle, feline, retail_milk, other |

## Processed data

### `data/processed/corpus_variants.parquet`

All concordant variants (iVar + LoFreq agreement) across the corpus. 289,054 rows.

| Variable | Type | Description |
|----------|------|-------------|
| `sample` | string | SRA run accession |
| `chrom` | string | Segment accession (e.g., PP755964.1) |
| `pos` | int | 1-indexed nucleotide position on segment |
| `ref` | string | Reference allele |
| `alt` | string | Alternative allele |
| `af_ivar` | float | Allele frequency reported by iVar |
| `af_lofreq` | float | Allele frequency reported by LoFreq |
| `af_mean` | float | Mean of iVar and LoFreq AF estimates |
| `depth_ivar` | int | Read depth at position (iVar) |
| `depth_lofreq` | int | Read depth at position (LoFreq) |
| `passes_strand_bias` | bool | True if Fisher's exact test p > 0.001 |

### `data/processed/corpus_variants_3pct.parquet`

Subset of `corpus_variants.parquet` with AF >= 3% and passing strand-bias filter. 214,226 rows.

### `data/processed/corpus_variants_subconsensus.parquet`

Subset of `corpus_variants_3pct.parquet` with AF between 3% and 50%. 56,933 rows.

### `data/processed/corpus_coverage.parquet`

Per-sample coverage summary. 4,559 rows.

| Variable | Type | Description |
|----------|------|-------------|
| `run_accession` | string | SRA run accession |
| `host_category` | string | cattle, feline, or retail_milk |
| `platform` | string | ILLUMINA or OXFORD_NANOPORE |
| `median_depth` | float | Median depth across all positions |
| `ivar_snvs` | int | Number of raw iVar variant calls |

### `data/processed/site_depth_by_sample.parquet`

Per-sample minimum codon depth at each of the 11 Tier 1 adaptation sites. 4,559 rows.

| Variable | Type | Description |
|----------|------|-------------|
| `sample` | string | SRA run accession |
| `S01_min_depth` | int | Min depth across PB2-627 codon positions (nt 1879-1881) |
| `S02_min_depth` | int | Min depth across PB2-701 codon positions (nt 2101-2103) |
| `S03_min_depth` | int | Min depth across PB2-591 codon positions (nt 1771-1773) |
| `S04_min_depth` | int | Min depth across PB2-271 codon positions (nt 811-813) |
| `S05_min_depth` | int | Min depth across PB2-631 codon positions (nt 1891-1893) |
| `S06_min_depth` | int | Min depth across PA-497 codon positions (nt 1489-1491) |
| `S07_min_depth` | int | Min depth across HA-238 codon positions (nt 712-714) |
| `S08_min_depth` | int | Min depth across HA-240 codon positions (nt 718-720) |
| `S09_min_depth` | int | Min depth across PB1-F2-66 codon positions (nt 290-292) |
| `S10_min_depth` | int | Min depth across NS1-92 codon positions (nt 274-276) |
| `S11_min_depth` | int | Min depth across M2-31 codon positions (nt 779-781) |

### `data/processed/segment_coverage_by_sample.parquet`

Per-sample, per-segment position counts at >= 100x depth. 4,559 rows.

| Variable | Type | Description |
|----------|------|-------------|
| `sample` | string | SRA run accession |
| `{SEGMENT}_covered` | int | Positions with >= 100x depth on that segment |
| `{SEGMENT}_total` | int | Total positions on that segment |
| `segments_with_coverage` | int | Number of segments with any position >= 100x |
| `genome_covered` | int | Total positions >= 100x across all segments |
| `genome_total` | int | Total positions across all segments (13,136) |

Where `{SEGMENT}` is one of: PB2, PB1, PA, HA, NP, NA, M, NS.

## Reference genome

**A/cattle/Texas/24-009308-004 (H5N1)**, B3.13 genotype. GenBank accessions:

| Segment | Accession | Length (nt) |
|---------|-----------|-------------|
| PB2 | PP755964.1 | 2,280 |
| PB1 | PP755963.1 | 2,274 |
| PA | PP755962.1 | 2,151 |
| HA | PP755957.1 | 1,704 |
| NP | PP755960.1 | 1,497 |
| NA | PP755959.1 | 1,410 |
| M | PP755958.1 | 982 |
| NS | PP755961.1 | 838 |
| **Total** | | **13,136** |
