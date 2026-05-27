#!/usr/bin/env bash
# MilkVar Phase 3: Full corpus processing on Vast.ai
# Processes all samples through the dual-caller iSNV pipeline using GNU parallel.
#
# Usage:
#   1. Upload manifest.tsv and this script to the Vast.ai instance
#   2. bash run_corpus_vastai.sh [--workers N] [--resume] manifest.tsv
#
# Output: results/ directory with variants/, coverage/, qc/, logs/
# Retrieve results: rsync -avz results/ local:path/to/results/
set -euo pipefail

WORKERS=12
RESUME=false
MANIFEST=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --workers) WORKERS="$2"; shift 2 ;;
        --resume) RESUME=true; shift ;;
        *) MANIFEST="$1"; shift ;;
    esac
done

if [[ -z "${MANIFEST}" || ! -f "${MANIFEST}" ]]; then
    echo "Usage: $0 [--workers N] [--resume] manifest.tsv"
    exit 1
fi

MANIFEST="$(cd "$(dirname "${MANIFEST}")" && pwd)/$(basename "${MANIFEST}")"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="/workspace/milkvar"
RESULTS="${WORK_DIR}/results"
REF="${WORK_DIR}/reference.fasta"
PROGRESS_FILE="${RESULTS}/completed.txt"
FAILED_FILE="${RESULTS}/failed.txt"

mkdir -p "${RESULTS}"/{variants,coverage,qc,logs}
touch "${PROGRESS_FILE}" "${FAILED_FILE}"

echo "============================================"
echo "MilkVar Corpus Processing — Vast.ai"
echo "============================================"
echo "Manifest:  ${MANIFEST}"
echo "Workers:   ${WORKERS}"
echo "Results:   ${RESULTS}"
echo "Resume:    ${RESUME}"
echo ""

# ── Install tools ──────────────────────────────────────────────────
install_tools() {
    echo "[1/4] Installing bioinformatics tools..."
    export DEBIAN_FRONTEND=noninteractive

    apt-get update -qq
    apt-get install -y -qq samtools parallel pigz python3-pip curl 2>/dev/null

    if ! command -v minimap2 &>/dev/null; then
        curl -sL https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
            | tar -xjf - -C /tmp
        cp /tmp/minimap2-2.28_x64-linux/minimap2 /usr/local/bin/
    fi

    if ! command -v fastp &>/dev/null; then
        curl -sL http://opengene.org/fastp/fastp.0.23.4 -o /usr/local/bin/fastp
        chmod +x /usr/local/bin/fastp
    fi

    if ! command -v ivar &>/dev/null; then
        apt-get install -y -qq autoconf automake libtool pkg-config \
            zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev 2>/dev/null
        (cd /tmp && curl -sL https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 | tar -xjf - \
            && cd htslib-1.20 && ./configure --prefix=/usr/local -q && make -j"$(nproc)" -s && make install -s && ldconfig)
        (cd /tmp && curl -sL https://github.com/andersen-lab/ivar/archive/refs/tags/v1.4.3.tar.gz | tar -xzf - \
            && cd ivar-1.4.3 && ./autogen.sh && ./configure --prefix=/usr/local --with-hts=/usr/local -q \
            && make -j"$(nproc)" -s && make install -s && ldconfig)
    fi

    if ! command -v lofreq &>/dev/null; then
        curl -sL https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz \
            | tar -xzf - -C /usr/local --strip-components=1
    fi

    if ! command -v prefetch &>/dev/null; then
        curl -sL https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
            | tar -xzf - -C /opt
        ln -sf /opt/sratoolkit.*/bin/prefetch /usr/local/bin/prefetch
        ln -sf /opt/sratoolkit.*/bin/fasterq-dump /usr/local/bin/fasterq-dump
        ln -sf /opt/sratoolkit.*/bin/vdb-config /usr/local/bin/vdb-config
    fi

    vdb-config --set /LIBS/GUID=vastai-milkvar \
        --set /libs/cloud/report_instance_identity=false 2>/dev/null || true

    pip install -q biopython pyarrow pandas tqdm 2>/dev/null

    echo "  Tools installed."
    for tool in minimap2 samtools ivar lofreq fastp prefetch fasterq-dump parallel; do
        printf "  %-20s %s\n" "${tool}:" "$(command -v ${tool} &>/dev/null && echo OK || echo MISSING)"
    done
}

# ── Download reference genome ──────────────────────────────────────
fetch_reference() {
    if [[ -f "${REF}" ]]; then
        echo "[2/4] Reference already exists."
        return
    fi
    echo "[2/4] Downloading reference genome..."
    python3 - <<'PYEOF'
from Bio import Entrez, SeqIO
Entrez.email = 'hayden.farquhar@icloud.com'
accs = ['PP755964.1','PP755963.1','PP755962.1','PP755957.1',
        'PP755960.1','PP755959.1','PP755958.1','PP755961.1']
handle = Entrez.efetch(db='nucleotide', id=','.join(accs), rettype='fasta', retmode='text')
with open('reference.fasta', 'w') as f:
    f.write(handle.read())
handle.close()
PYEOF
    minimap2 -d "${REF}.mmi" "${REF}" 2>/dev/null
    echo "  Reference ready: $(grep -c '^>' "${REF}") segments"
}

# ── Build sample list ──────────────────────────────────────────────
build_worklist() {
    echo "[3/4] Building worklist..."
    TOTAL=$(tail -n +2 "${MANIFEST}" | wc -l)

    if ${RESUME} && [[ -s "${PROGRESS_FILE}" ]]; then
        DONE=$(wc -l < "${PROGRESS_FILE}")
        echo "  Resuming: ${DONE}/${TOTAL} already completed"
    else
        > "${PROGRESS_FILE}"
        > "${FAILED_FILE}"
    fi

    WORKLIST="${WORK_DIR}/worklist.tsv"
    tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r acc rest; do
        if ! grep -qx "${acc}" "${PROGRESS_FILE}" 2>/dev/null; then
            echo "${acc}"
        fi
    done > "${WORK_DIR}/remaining_accs.txt"

    # Build worklist with platform info
    python3 - "${MANIFEST}" "${WORK_DIR}/remaining_accs.txt" "${WORKLIST}" <<'PYEOF'
import sys, csv
manifest_path, remaining_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
remaining = set(open(remaining_path).read().strip().split('\n'))
with open(manifest_path) as f, open(out_path, 'w') as out:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        acc = row['run_accession']
        if acc in remaining:
            platform = row.get('platform', 'ILLUMINA')
            library = row.get('library_strategy', 'WGS')
            out.write(f"{acc}\t{platform}\t{library}\n")
PYEOF

    REMAINING=$(wc -l < "${WORKLIST}")
    echo "  Worklist: ${REMAINING} samples to process"
}

# ── Process wrapper (called by parallel) ───────────────────────────
export -f 2>/dev/null || true
export REF RESULTS PROGRESS_FILE FAILED_FILE

process_one() {
    local acc="$1" platform="$2" library="$3"

    if bash "${SCRIPT_DIR}/process_sample.sh" "${acc}" "${platform}" "${library}" "${REF}" "${RESULTS}"; then
        echo "${acc}" >> "${PROGRESS_FILE}"
    else
        echo "${acc}" >> "${FAILED_FILE}"
        echo "FAIL: ${acc}" >&2
    fi
}
export -f process_one
export SCRIPT_DIR

# ── Main processing loop ──────────────────────────────────────────
run_parallel() {
    echo "[4/4] Processing with ${WORKERS} parallel workers..."
    echo ""

    WORKLIST="${WORK_DIR}/worklist.tsv"
    TOTAL_REMAINING=$(wc -l < "${WORKLIST}")

    START_TIME=$(date +%s)

    cat "${WORKLIST}" | parallel --colsep '\t' -j "${WORKERS}" --halt soon,fail=20% \
        process_one {1} {2} {3}

    END_TIME=$(date +%s)
    ELAPSED=$(( END_TIME - START_TIME ))
    COMPLETED=$(wc -l < "${PROGRESS_FILE}")
    FAILED=$(wc -l < "${FAILED_FILE}")

    echo ""
    echo "============================================"
    echo "PROCESSING COMPLETE"
    echo "============================================"
    echo "  Completed:  ${COMPLETED}"
    echo "  Failed:     ${FAILED}"
    echo "  Wall time:  $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"
    echo "  Rate:       $(python3 -c "print(f'{${TOTAL_REMAINING}/(${ELAPSED}/60):.1f} samples/min')" 2>/dev/null || echo "N/A")"
    echo ""
    echo "Results in: ${RESULTS}/"
    echo "  variants/  — *.ivar.tsv + *.lofreq.vcf"
    echo "  coverage/  — *.depth.tsv.gz"
    echo "  qc/        — *.fastp.json"
    echo "  logs/      — per-sample logs"
    echo ""
    echo "Next: rsync results back and run concordance_filter.py"
}

# ── Entry point ────────────────────────────────────────────────────
install_tools
fetch_reference
build_worklist
run_parallel
