#!/usr/bin/env bash
# MilkVar Phase 4: C4 Positive-Control Validation
# V1: Consensus recovery vs GenBank (5 smoke-test samples)
# V3: Synthetic spike-in at 4 mixture ratios
#
# Usage: bash run_validation.sh
set -euo pipefail

WORK="/workspace/milkvar_validation"
REF="${WORK}/reference.fasta"
mkdir -p "${WORK}"/{v1_consensus,v3_spikein,bams,logs}

# ── Install tools ──────────────────────────────────────────────────
echo "[1/5] Installing tools..."
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
    vdb-config --set /LIBS/GUID=vastai-validation \
        --set /libs/cloud/report_instance_identity=false 2>/dev/null || true
fi

pip install -q biopython pyarrow pandas 2>/dev/null
echo "  Tools ready."

# ── Download reference ─────────────────────────────────────────────
echo "[2/5] Downloading reference..."
if [[ ! -f "${REF}" ]]; then
    python3 -c "
from Bio import Entrez
Entrez.email = 'hayden.farquhar@icloud.com'
accs = ['PP755964.1','PP755963.1','PP755962.1','PP755957.1',
        'PP755960.1','PP755959.1','PP755958.1','PP755961.1']
h = Entrez.efetch(db='nucleotide', id=','.join(accs), rettype='fasta', retmode='text')
open('${REF}', 'w').write(h.read())
h.close()
"
    minimap2 -d "${REF}.mmi" "${REF}" 2>/dev/null
fi
echo "  Reference ready."

# ── Helper: download + align one sample ────────────────────────────
fetch_and_align() {
    local acc="$1"
    local prefix="${acc:0:6}"
    local sub=""
    [[ ${#acc} -gt 9 ]] && sub="/$(printf "%03d" "${acc:9}")"
    local base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}${sub}/${acc}"
    local tmpdir="/tmp/val_${acc}"
    mkdir -p "${tmpdir}"

    # Download from ENA
    for suffix in _1.fastq.gz _2.fastq.gz .fastq.gz; do
        curl -sS -f --retry 2 --connect-timeout 15 --max-time 300 \
            -o "${tmpdir}/${acc}${suffix}" "${base}/${acc}${suffix}" 2>/dev/null || rm -f "${tmpdir}/${acc}${suffix}"
    done

    # SRA fallback
    if [[ ! -f "${tmpdir}/${acc}_1.fastq.gz" && ! -f "${tmpdir}/${acc}.fastq.gz" ]]; then
        prefetch "${acc}" -O "${tmpdir}/" --max-size 50G 2>/dev/null || true
        local sra="${tmpdir}/${acc}/${acc}.sra"
        [[ ! -f "${sra}" ]] && sra="${tmpdir}/${acc}.sra"
        [[ -f "${sra}" ]] && fasterq-dump "${sra}" -O "${tmpdir}/" --split-files --threads 4 2>/dev/null && \
            gzip -1 "${tmpdir}"/${acc}*.fastq
    fi

    # Align
    local bam="${WORK}/bams/${acc}.sorted.bam"
    if [[ -f "${tmpdir}/${acc}_1.fastq.gz" ]]; then
        minimap2 -a -x sr -t 4 "${REF}" "${tmpdir}/${acc}_1.fastq.gz" "${tmpdir}/${acc}_2.fastq.gz" 2>/dev/null \
            | samtools view -bS -q 20 - | samtools sort -@ 4 -o "${bam}"
    elif [[ -f "${tmpdir}/${acc}.fastq.gz" ]]; then
        minimap2 -a -x sr -t 4 "${REF}" "${tmpdir}/${acc}.fastq.gz" 2>/dev/null \
            | samtools view -bS -q 20 - | samtools sort -@ 4 -o "${bam}"
    fi
    samtools index "${bam}"
    rm -rf "${tmpdir}"
    echo "${bam}"
}

# ── V1: Consensus recovery ────────────────────────────────────────
echo "[3/5] V1: Consensus recovery for 5 samples..."

V1_SAMPLES=("SRR35150577" "SRR29165098" "SRR36114813" "SRR38788355" "SRR29165088")
V1_DIR="${WORK}/v1_consensus"

for acc in "${V1_SAMPLES[@]}"; do
    echo "  Downloading & aligning ${acc}..."
    bam=$(fetch_and_align "${acc}")

    # Extract consensus
    samtools mpileup -aa -A -d 8000 -Q 0 --reference "${REF}" "${bam}" 2>/dev/null \
        | ivar consensus -p "${V1_DIR}/${acc}.consensus" -q 20 -t 0.5 -m 20 2>/dev/null

    echo "  ${acc}: consensus extracted"
done

# Download GenBank consensus for comparison
echo "  Downloading GenBank consensus sequences for comparison..."
python3 - "${V1_DIR}" <<'PYEOF'
import sys, os
from Bio import Entrez, SeqIO
Entrez.email = 'hayden.farquhar@icloud.com'

v1_dir = sys.argv[1]
samples = ['SRR35150577', 'SRR29165098', 'SRR36114813', 'SRR38788355', 'SRR29165088']

# For each sample, try to find the linked GenBank consensus via BioSample
for acc in samples:
    try:
        handle = Entrez.esearch(db='sra', term=acc)
        result = Entrez.read(handle)
        handle.close()
        if result['IdList']:
            handle = Entrez.elink(dbfrom='sra', db='nuccore', id=result['IdList'][0])
            links = Entrez.read(handle)
            handle.close()
            nuccore_ids = []
            for linkset in links:
                for link in linkset.get('LinkSetDb', []):
                    nuccore_ids.extend([l['Id'] for l in link['Link']])
            if nuccore_ids:
                handle = Entrez.efetch(db='nuccore', id=nuccore_ids[:8], rettype='fasta', retmode='text')
                with open(os.path.join(v1_dir, f'{acc}.genbank.fasta'), 'w') as f:
                    f.write(handle.read())
                handle.close()
                print(f'  {acc}: GenBank consensus downloaded ({len(nuccore_ids)} segments)')
            else:
                print(f'  {acc}: No linked GenBank sequences found')
        else:
            print(f'  {acc}: SRA record not found')
    except Exception as e:
        print(f'  {acc}: Error - {e}')
PYEOF

# Compare consensus
echo "  Comparing consensus..."
python3 - "${V1_DIR}" "${REF}" <<'PYEOF'
import sys, os
from Bio import SeqIO

v1_dir, ref_path = sys.argv[1], sys.argv[2]
samples = ['SRR35150577', 'SRR29165098', 'SRR36114813', 'SRR38788355', 'SRR29165088']

ref_seqs = {r.id: str(r.seq) for r in SeqIO.parse(ref_path, 'fasta')}
print(f'\n=== V1 CONSENSUS RECOVERY ===')
print(f'Reference segments: {len(ref_seqs)}')

for acc in samples:
    cons_path = os.path.join(v1_dir, f'{acc}.consensus.fa')
    if not os.path.exists(cons_path):
        print(f'{acc}: NO CONSENSUS FILE')
        continue

    cons_seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(cons_path, 'fasta')}

    total_pos = 0
    mismatches = 0
    n_positions = 0
    for seg_id, ref_seq in ref_seqs.items():
        cons_seq = cons_seqs.get(seg_id, '')
        if not cons_seq:
            continue
        min_len = min(len(ref_seq), len(cons_seq))
        for i in range(min_len):
            r, c = ref_seq[i], cons_seq[i]
            if c == 'N' or c == 'n':
                continue
            n_positions += 1
            if r != c:
                mismatches += 1

    print(f'{acc}: {mismatches} mismatches across {n_positions} positions')
print()
PYEOF

# ── V3: Synthetic spike-in ────────────────────────────────────────
echo "[4/5] V3: Synthetic spike-in..."

# Use two high-coverage cattle samples from different states
SAMPLE_A="SRR29165118"  # Texas cattle, ~7420x
SAMPLE_B="SRR29165114"  # Ohio cattle, ~5277x

echo "  Downloading ${SAMPLE_A} (Texas) and ${SAMPLE_B} (Ohio)..."
BAM_A=$(fetch_and_align "${SAMPLE_A}")
BAM_B=$(fetch_and_align "${SAMPLE_B}")

echo "  Finding consensus-level differences..."
V3_DIR="${WORK}/v3_spikein"

python3 - "${BAM_A}" "${BAM_B}" "${REF}" "${V3_DIR}" <<'PYEOF'
import subprocess, sys, os
from collections import defaultdict

bam_a, bam_b, ref, outdir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
os.makedirs(outdir, exist_ok=True)

def get_consensus(bam, ref):
    """Extract per-position consensus from BAM."""
    consensus = {}
    cmd = f'samtools mpileup -aa -A -d 8000 -Q 20 --reference {ref} {bam}'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    for line in proc.stdout.split('\n'):
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        chrom, pos, ref_base, depth = parts[0], int(parts[1]), parts[2], int(parts[3])
        if depth < 50:
            continue
        pileup = parts[4].upper()
        counts = defaultdict(int)
        counts[ref_base.upper()] = pileup.count('.') + pileup.count(',')
        for base in 'ACGT':
            if base != ref_base.upper():
                counts[base] = pileup.count(base)
        total = sum(counts.values())
        if total == 0:
            continue
        major = max(counts, key=counts.get)
        major_freq = counts[major] / total
        if major_freq >= 0.75:
            consensus[(chrom, pos)] = (major, major_freq, depth)
    return consensus

print('  Extracting consensus for sample A...')
cons_a = get_consensus(bam_a, ref)
print(f'    {len(cons_a)} positions')

print('  Extracting consensus for sample B...')
cons_b = get_consensus(bam_b, ref)
print(f'    {len(cons_b)} positions')

# Find positions where A and B differ at consensus level
diffs = []
shared_positions = set(cons_a.keys()) & set(cons_b.keys())
for pos in sorted(shared_positions):
    base_a, freq_a, dp_a = cons_a[pos]
    base_b, freq_b, dp_b = cons_b[pos]
    if base_a != base_b and freq_a >= 0.90 and freq_b >= 0.90:
        diffs.append({
            'chrom': pos[0], 'pos': pos[1],
            'base_a': base_a, 'freq_a': freq_a, 'depth_a': dp_a,
            'base_b': base_b, 'freq_b': freq_b, 'depth_b': dp_b,
        })

print(f'\n  Consensus-level differences (both >=90% AF): {len(diffs)}')
if diffs:
    for d in diffs[:20]:
        print(f"    {d['chrom']}:{d['pos']}  A={d['base_a']}({d['freq_a']:.1%})  B={d['base_b']}({d['freq_b']:.1%})")

# Save truth set
with open(os.path.join(outdir, 'truth_set.tsv'), 'w') as f:
    f.write('chrom\tpos\tbase_a\tfreq_a\tbase_b\tfreq_b\n')
    for d in diffs:
        f.write(f"{d['chrom']}\t{d['pos']}\t{d['base_a']}\t{d['freq_a']:.4f}\t{d['base_b']}\t{d['freq_b']:.4f}\n")

print(f'\n  Truth set saved: {len(diffs)} positions')
if len(diffs) < 5:
    print('  WARNING: Need >= 5 differences for spike-in. Trying different samples may be needed.')
PYEOF

# Create mixtures and run dual-caller pipeline
echo "  Creating read mixtures and running pipeline..."
python3 - "${BAM_A}" "${BAM_B}" "${REF}" "${V3_DIR}" <<'PYEOF'
import subprocess, os, sys
import pandas as pd

bam_a, bam_b, ref, outdir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# Load truth set
truth = pd.read_csv(os.path.join(outdir, 'truth_set.tsv'), sep='\t')
if len(truth) < 5:
    print('  Insufficient differences for spike-in. Exiting.')
    sys.exit(0)

truth_keys = set(truth['chrom'] + ':' + truth['pos'].astype(str) + ':' + truth['base_b'])
print(f'  Truth set: {len(truth_keys)} spiked-in variant positions')

# Mixture ratios: A is major, B is minor (spiked in)
ratios = [(0.95, 0.05, '95_05'), (0.99, 0.01, '99_01'),
          (0.995, 0.005, '995_005'), (0.998, 0.002, '998_002')]

# Get total read counts
count_a = int(subprocess.run(f'samtools view -c {bam_a}', shell=True, capture_output=True, text=True).stdout.strip())
count_b = int(subprocess.run(f'samtools view -c {bam_b}', shell=True, capture_output=True, text=True).stdout.strip())
target_total = min(count_a, count_b, 200000)

print(f'  Reads: A={count_a}, B={count_b}, target_total={target_total}')

results = []
for frac_a, frac_b, label in ratios:
    print(f'\n  Mixture {label} (A:{frac_a*100:.1f}% / B:{frac_b*100:.1f}%)...')
    mix_bam = os.path.join(outdir, f'mix_{label}.sorted.bam')

    n_a = int(target_total * frac_a)
    n_b = int(target_total * frac_b)
    seed_a = 42
    seed_b = 43
    subsample_a = n_a / count_a if n_a < count_a else 1.0
    subsample_b = n_b / count_b if n_b < count_b else 1.0

    # Subsample and merge
    tmp_a = os.path.join(outdir, f'tmp_a_{label}.bam')
    tmp_b = os.path.join(outdir, f'tmp_b_{label}.bam')
    subprocess.run(f'samtools view -bs {seed_a}.{subsample_a:.6f} {bam_a} -o {tmp_a}', shell=True)
    subprocess.run(f'samtools view -bs {seed_b}.{subsample_b:.6f} {bam_b} -o {tmp_b}', shell=True)
    subprocess.run(f'samtools merge -f {mix_bam} {tmp_a} {tmp_b}', shell=True)
    subprocess.run(f'samtools sort -@ 4 -o {mix_bam}.tmp {mix_bam} && mv {mix_bam}.tmp {mix_bam}', shell=True)
    subprocess.run(f'samtools index {mix_bam}', shell=True)
    os.remove(tmp_a)
    os.remove(tmp_b)

    actual_count = int(subprocess.run(f'samtools view -c {mix_bam}', shell=True, capture_output=True, text=True).stdout.strip())
    print(f'    Mixed BAM: {actual_count} reads')

    # Run iVar
    ivar_prefix = os.path.join(outdir, f'mix_{label}.ivar')
    subprocess.run(
        f'samtools mpileup -aa -A -d 8000 -B -Q 0 --reference {ref} {mix_bam} '
        f'| ivar variants -p {ivar_prefix} -q 20 -t 0.01 -m 10 -r {ref}',
        shell=True, capture_output=True)

    # Run LoFreq
    lofreq_out = os.path.join(outdir, f'mix_{label}.lofreq.vcf')
    subprocess.run(
        f'lofreq call --no-default-filter -f {ref} -o {lofreq_out} {mix_bam}',
        shell=True, capture_output=True, timeout=600)

    # Parse iVar results
    ivar_path = f'{ivar_prefix}.tsv'
    ivar_keys = set()
    if os.path.exists(ivar_path):
        try:
            idf = pd.read_csv(ivar_path, sep='\t')
            snvs = idf[(idf['REF'].str.len()==1) & (idf['ALT'].str.len()==1) & (~idf['ALT'].isin(['+','-']))]
            ivar_keys = set(snvs['REGION'] + ':' + snvs['POS'].astype(str) + ':' + snvs['ALT'])
        except: pass

    # Parse LoFreq results
    lofreq_keys = set()
    if os.path.exists(lofreq_out):
        with open(lofreq_out) as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) >= 5 and len(parts[3])==1 and len(parts[4])==1:
                    lofreq_keys.add(f'{parts[0]}:{parts[1]}:{parts[4]}')

    concordant = ivar_keys & lofreq_keys

    # Compute sensitivity and specificity at each threshold
    for thresh_label, thresh_val in [('1%', 0.01), ('2%', 0.02), ('3%', 0.03), ('5%', 0.05)]:
        # Filter iVar by threshold
        ivar_at_t = set()
        if os.path.exists(ivar_path):
            try:
                idf2 = pd.read_csv(ivar_path, sep='\t')
                snvs2 = idf2[(idf2['ALT_FREQ']>=thresh_val) & (idf2['REF'].str.len()==1) & (idf2['ALT'].str.len()==1) & (~idf2['ALT'].isin(['+','-']))]
                ivar_at_t = set(snvs2['REGION'] + ':' + snvs2['POS'].astype(str) + ':' + snvs2['ALT'])
            except: pass

        lofreq_at_t = set()
        if os.path.exists(lofreq_out):
            with open(lofreq_out) as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 8 and len(parts[3])==1 and len(parts[4])==1:
                        info = dict(x.split('=',1) for x in parts[7].split(';') if '=' in x)
                        if float(info.get('AF',0)) >= thresh_val:
                            lofreq_at_t.add(f'{parts[0]}:{parts[1]}:{parts[4]}')

        conc_at_t = ivar_at_t & lofreq_at_t

        tp = len(conc_at_t & truth_keys)
        fn = len(truth_keys - conc_at_t)
        fp = len(conc_at_t - truth_keys)

        # Total non-spiked positions for specificity
        # Approximate: genome is 13136 bp, truth has N positions
        total_neg = 13136 - len(truth_keys)
        sensitivity = tp / len(truth_keys) if truth_keys else 0
        specificity = 1 - (fp / total_neg) if total_neg > 0 else 0

        results.append({
            'mixture': label, 'expected_af': frac_b,
            'threshold': thresh_label, 'thresh_val': thresh_val,
            'tp': tp, 'fn': fn, 'fp': fp,
            'sensitivity': sensitivity, 'specificity': specificity,
            'ivar_total': len(ivar_at_t), 'lofreq_total': len(lofreq_at_t),
            'concordant_total': len(conc_at_t),
        })

# Print results
rdf = pd.DataFrame(results)
rdf.to_csv(os.path.join(outdir, 'spikein_results.tsv'), sep='\t', index=False)

print(f'\n{"="*90}')
print(f'{"V3 SYNTHETIC SPIKE-IN RESULTS":^90}')
print(f'{"="*90}')
print(f'Truth set: {len(truth_keys)} consensus-level SNV differences between samples A and B')
print()
print(f'{"Mixture":>10} {"Exp AF":>8} {"Thresh":>8} {"TP":>5} {"FN":>5} {"FP":>5} {"Sensitivity":>12} {"Specificity":>12} {"PASS?":>7}')
print('-' * 80)
for _, r in rdf.iterrows():
    passed = 'YES' if r['sensitivity'] >= 0.80 and r['specificity'] >= 0.99 else ''
    if r['sensitivity'] >= 0.80:
        passed = 'SENS' if r['specificity'] < 0.99 else 'YES'
    print(f'{r["mixture"]:>10} {r["expected_af"]:>7.1%} {r["threshold"]:>8} {r["tp"]:>5} {r["fn"]:>5} {r["fp"]:>5} {r["sensitivity"]:>11.1%} {r["specificity"]:>11.4%} {passed:>7}')
PYEOF

# ── Summary ────────────────────────────────────────────────────────
echo ""
echo "[5/5] Validation complete."
echo "Results in: ${WORK}/"
echo "  v1_consensus/ — per-sample consensus FASTA + GenBank comparisons"
echo "  v3_spikein/   — mixture BAMs, variant calls, truth set, results table"
