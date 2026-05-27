#!/usr/bin/env bash
# Process a single H5N1 sample through: fetch → QC → align → depth → iVar → LoFreq
# Called by GNU parallel from run_corpus_vastai.sh
# Usage: process_sample.sh <accession> <platform> <library> <ref_fasta> <results_dir>
set -euo pipefail

ACC="$1"
PLATFORM="$2"
LIBRARY="$3"
REF="$4"
RESULTS="$5"

VARIANTS="${RESULTS}/variants"
COVERAGE="${RESULTS}/coverage"
QC="${RESULTS}/qc"
TMPDIR="/tmp/milkvar/${ACC}"
LOGFILE="${RESULTS}/logs/${ACC}.log"

mkdir -p "${TMPDIR}" "${VARIANTS}" "${COVERAGE}" "${QC}" "${RESULTS}/logs"

exec > "${LOGFILE}" 2>&1

cleanup() { rm -rf "${TMPDIR}"; }
trap cleanup EXIT

MIN_MAPQ=20
MIN_FREQ=0.01
MIN_DEPTH=10
DEPTH_CAP=8000
DOWNSAMPLE_THRESH=10000

is_nanopore=false
[[ "${PLATFORM}" == "OXFORD_NANOPORE" ]] && is_nanopore=true

# ── Stage 1: Fetch FASTQs (ENA first, SRA fallback) ───────────────
ena_download() {
    local acc="$1" outdir="$2"
    local prefix="${acc:0:6}"
    local sub=""
    if [[ ${#acc} -gt 9 ]]; then
        sub="/$(printf "%03d" "${acc:9}")"
    fi
    local base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}${sub}/${acc}"

    # Try paired first, then single
    local got_r1=false got_r2=false got_single=false
    for suffix in _1.fastq.gz _2.fastq.gz .fastq.gz; do
        local url="${base}/${acc}${suffix}"
        local out="${outdir}/${acc}${suffix}"
        if curl -sS -f --retry 2 --retry-delay 5 --connect-timeout 15 --max-time 300 \
                -o "${out}" "${url}" 2>/dev/null; then
            case "${suffix}" in
                _1.fastq.gz) got_r1=true ;;
                _2.fastq.gz) got_r2=true ;;
                .fastq.gz)   got_single=true ;;
            esac
        else
            rm -f "${out}"
        fi
    done
    ${got_r1} || ${got_single}
}

sra_download() {
    local acc="$1" outdir="$2"
    prefetch "${acc}" -O "${outdir}/" --max-size 50G 2>/dev/null || true
    local sra_file="${outdir}/${acc}/${acc}.sra"
    [[ ! -f "${sra_file}" ]] && sra_file="${outdir}/${acc}.sra"
    if [[ -f "${sra_file}" ]]; then
        fasterq-dump "${sra_file}" -O "${outdir}/" --split-files --threads 2 2>/dev/null || true
        rm -f "${sra_file}"
        rmdir "${outdir}/${acc}" 2>/dev/null || true
        # Compress to match ENA format so downstream is uniform
        for fq in "${outdir}"/${acc}*.fastq; do
            [[ -f "${fq}" ]] && gzip -1 "${fq}" &
        done
        wait
    fi
}

echo "DOWNLOAD: trying ENA..."
if ! ena_download "${ACC}" "${TMPDIR}"; then
    echo "DOWNLOAD: ENA failed, trying SRA..."
    sra_download "${ACC}" "${TMPDIR}"
fi

# Identify files (all .gz at this point)
R1_GZ="${TMPDIR}/${ACC}_1.fastq.gz"
R2_GZ="${TMPDIR}/${ACC}_2.fastq.gz"
SINGLE_GZ="${TMPDIR}/${ACC}.fastq.gz"

IS_PAIRED=false
if [[ -f "${R1_GZ}" && -f "${R2_GZ}" ]]; then
    IS_PAIRED=true
elif [[ ! -f "${SINGLE_GZ}" ]]; then
    echo "FAIL: no FASTQ data" >&2
    exit 1
fi

# ── Stage 2: fastp QC (reads .gz directly) ─────────────────────────
if ${is_nanopore}; then QUAL=7; MINLEN=200; else QUAL=15; MINLEN=50; fi

if ${IS_PAIRED}; then
    fastp -i "${R1_GZ}" -I "${R2_GZ}" \
        -o "${TMPDIR}/${ACC}_qc_1.fastq" -O "${TMPDIR}/${ACC}_qc_2.fastq" \
        --qualified_quality_phred ${QUAL} --length_required ${MINLEN} \
        --detect_adapter_for_pe --thread 2 \
        --json "${QC}/${ACC}.fastp.json" --html /dev/null 2>/dev/null
    READS="${TMPDIR}/${ACC}_qc_1.fastq ${TMPDIR}/${ACC}_qc_2.fastq"
else
    fastp -i "${SINGLE_GZ}" \
        -o "${TMPDIR}/${ACC}_qc.fastq" \
        --qualified_quality_phred ${QUAL} --length_required ${MINLEN} \
        --thread 2 \
        --json "${QC}/${ACC}.fastp.json" --html /dev/null 2>/dev/null
    READS="${TMPDIR}/${ACC}_qc.fastq"
fi

# Free raw FASTQs immediately
rm -f "${R1_GZ}" "${R2_GZ}" "${SINGLE_GZ}" 2>/dev/null

# ── Stage 3: Alignment ─────────────────────────────────────────────
if ${is_nanopore}; then PRESET="map-ont"; else PRESET="sr"; fi
BAM="${TMPDIR}/${ACC}.sorted.bam"

minimap2 -a -x "${PRESET}" -t 2 "${REF}" ${READS} 2>/dev/null \
    | samtools view -bS -q ${MIN_MAPQ} - \
    | samtools sort -@ 2 -o "${BAM}"
samtools index "${BAM}"

# Free QC'd FASTQs
rm -f ${READS} 2>/dev/null

# ── Stage 3b: Downsample if extremely deep ─────────────────────────
GENOME_LEN=13136
MAPPED=$(samtools idxstats "${BAM}" | awk '{s+=$3} END {print s}')
READ_LEN=$(samtools stats "${BAM}" 2>/dev/null | awk '/^SN\taverage length:/ {print int($NF)}')
[[ -z "${READ_LEN}" || "${READ_LEN}" -eq 0 ]] && READ_LEN=150
EST_DEPTH=$(( MAPPED * READ_LEN / GENOME_LEN ))

if [[ "${EST_DEPTH}" -gt "${DOWNSAMPLE_THRESH}" ]]; then
    FRAC=$(python3 -c "print(f'{${DOWNSAMPLE_THRESH}/${EST_DEPTH}:.4f}')")
    samtools view -bs "42.${FRAC#0.}" "${BAM}" | samtools sort -@ 2 -o "${TMPDIR}/${ACC}.ds.bam"
    mv "${TMPDIR}/${ACC}.ds.bam" "${BAM}"
    samtools index "${BAM}"
    echo "DOWNSAMPLED: ${EST_DEPTH}x → ~${DOWNSAMPLE_THRESH}x (frac=${FRAC})"
fi

# ── Stage 4: Depth ──────────────────────────────────────────────────
samtools depth -a "${BAM}" | gzip > "${COVERAGE}/${ACC}.depth.tsv.gz"

# ── Stage 5: iVar ──────────────────────────────────────────────────
timeout 600 bash -c "samtools mpileup -aa -A -d ${DEPTH_CAP} -B -Q 0 --reference ${REF} ${BAM} \
    | ivar variants -p ${VARIANTS}/${ACC}.ivar -q 20 -t ${MIN_FREQ} -m ${MIN_DEPTH} -r ${REF}"

# ── Stage 6: LoFreq ────────────────────────────────────────────────
if ${is_nanopore}; then
    printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
        > "${VARIANTS}/${ACC}.lofreq.vcf"
else
    timeout 900 lofreq call --no-default-filter -f "${REF}" \
        -o "${VARIANTS}/${ACC}.lofreq.vcf" "${BAM}" 2>/dev/null || \
        printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
            > "${VARIANTS}/${ACC}.lofreq.vcf"
fi

echo "OK"
