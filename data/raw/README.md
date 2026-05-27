# Raw Data Acquisition

The raw sequencing reads are publicly available on the NCBI Sequence Read Archive. They are not included in this repository due to size (>500 GB total).

## Manifest

`manifest.tsv` contains the complete corpus of 4,571 SRA accessions across 7 BioProjects. To download the reads:

```bash
# Install sra-tools
# https://github.com/ncbi/sra-tools

# Download a single sample
prefetch SRR28752447
fasterq-dump SRR28752447

# Or use ENA for faster downloads (recommended)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR287/047/SRR28752447/SRR28752447_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR287/047/SRR28752447/SRR28752447_2.fastq.gz
```

## Reference genome

Download the B3.13 reference genome:

```bash
# Fetch all 8 segments
for acc in PP755957 PP755958 PP755959 PP755960 PP755961 PP755962 PP755963 PP755964; do
    efetch -db nucleotide -id ${acc}.1 -format fasta >> reference.fasta
done
```

## Pipeline dependencies

The full variant-calling pipeline (scripts 02-07) requires:

- BWA (0.7.18+)
- minimap2 (2.28+)
- samtools (1.21+)
- iVar (1.4.3+)
- LoFreq (2.1.5+)
- fastp (0.23.4+)
- sra-tools (3.1.1+)
- GNU parallel

Install via conda:

```bash
conda create -n milkvar -c bioconda -c conda-forge \
    bwa minimap2 samtools ivar lofreq fastp sra-tools parallel
conda activate milkvar
```
