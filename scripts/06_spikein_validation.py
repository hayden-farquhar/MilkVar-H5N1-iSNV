#!/usr/bin/env python3
"""V3 Synthetic spike-in validation for MilkVar Phase 4."""
import subprocess, os, sys
import pandas as pd
from collections import defaultdict

REF = "validation/reference.fasta"
BAM_A = "validation/bams/spk_A.sorted.bam"
BAM_B = "validation/bams/spk_B.sorted.bam"
V3DIR = "validation/v3_spikein"
os.makedirs(V3DIR, exist_ok=True)

def get_consensus_bases(bam):
    consensus = {}
    cmd = f"samtools mpileup -aa -A -d 8000 -Q 20 --reference {REF} {bam}"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    for line in proc.stdout.split("\n"):
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, ref_base, depth = parts[0], int(parts[1]), parts[2].upper(), int(parts[3])
        if depth < 100:
            continue
        pileup = parts[4]
        counts = defaultdict(int)
        counts[ref_base] = pileup.count(".") + pileup.count(",")
        for b in "ACGT":
            if b != ref_base:
                counts[b] = pileup.upper().count(b)
        total = sum(counts.values())
        if total == 0:
            continue
        major = max(counts, key=counts.get)
        major_freq = counts[major] / total
        consensus[(chrom, pos)] = (major, major_freq, depth)
    return consensus

print("Finding consensus-level SNP differences...")
cons_a = get_consensus_bases(BAM_A)
cons_b = get_consensus_bases(BAM_B)
print(f"  A: {len(cons_a)} positions  |  B: {len(cons_b)} positions")

shared = set(cons_a.keys()) & set(cons_b.keys())
diffs = []
for pos in sorted(shared):
    ba, fa, da = cons_a[pos]
    bb, fb, db = cons_b[pos]
    if ba != bb and fa >= 0.90 and fb >= 0.90:
        diffs.append({"chrom": pos[0], "pos": pos[1],
                       "base_a": ba, "freq_a": fa,
                       "base_b": bb, "freq_b": fb})

print(f"  Consensus-level SNP differences: {len(diffs)}")
for d in diffs:
    print(f"    {d['chrom']}:{d['pos']}  A={d['base_a']}({d['freq_a']:.1%})  B={d['base_b']}({d['freq_b']:.1%})")

if len(diffs) < 5:
    print(f"\nWARNING: Only {len(diffs)} differences (need >=5).")
    if len(diffs) == 0:
        sys.exit(1)

truth_keys = set()
with open(os.path.join(V3DIR, "truth_set.tsv"), "w") as f:
    f.write("chrom\tpos\tbase_a\tbase_b\n")
    for d in diffs:
        f.write(f"{d['chrom']}\t{d['pos']}\t{d['base_a']}\t{d['base_b']}\n")
        truth_keys.add(f"{d['chrom']}:{d['pos']}:{d['base_b']}")

print(f"\nTruth set: {len(truth_keys)} positions where B differs from A")

count_a = int(subprocess.run(f"samtools view -c {BAM_A}", shell=True, capture_output=True, text=True).stdout.strip())
count_b = int(subprocess.run(f"samtools view -c {BAM_B}", shell=True, capture_output=True, text=True).stdout.strip())
target = min(count_a, count_b, 500000)

ratios = [(0.95, 0.05, "95_05"), (0.99, 0.01, "99_01"),
          (0.995, 0.005, "995_005"), (0.998, 0.002, "998_002")]

results = []
for frac_a, frac_b, label in ratios:
    print(f"\n--- Mixture {label} (A={frac_a*100:.1f}% / B={frac_b*100:.1f}%) ---")
    n_a = int(target * frac_a)
    n_b = int(target * frac_b)
    sub_a = min(n_a / count_a, 0.9999)
    sub_b = min(n_b / count_b, 0.9999)

    # samtools -s format: SEED.FRAC where FRAC is digits after decimal (no leading 0.)
    frac_str_a = f"{sub_a:.6f}"[2:]  # "0.298800" -> "298800"
    frac_str_b = f"{sub_b:.6f}"[2:]

    mix_bam = os.path.join(V3DIR, f"mix_{label}.bam")
    tmp_a = os.path.join(V3DIR, "tmp_a.bam")
    tmp_b = os.path.join(V3DIR, "tmp_b.bam")

    subprocess.run(f"samtools view -bs 42.{frac_str_a} {BAM_A} -o {tmp_a}", shell=True)
    subprocess.run(f"samtools view -bs 43.{frac_str_b} {BAM_B} -o {tmp_b}", shell=True)
    subprocess.run(f"samtools merge -f {mix_bam} {tmp_a} {tmp_b}", shell=True)
    sorted_bam = mix_bam.replace(".bam", ".sorted.bam")
    subprocess.run(f"samtools sort -@ 4 -o {sorted_bam} {mix_bam}", shell=True)
    subprocess.run(f"samtools index {sorted_bam}", shell=True)
    for f in [tmp_a, tmp_b, mix_bam]:
        if os.path.exists(f):
            os.remove(f)

    actual = int(subprocess.run(f"samtools view -c {sorted_bam}", shell=True, capture_output=True, text=True).stdout.strip())
    print(f"  Reads: {actual}")

    ivar_p = os.path.join(V3DIR, f"mix_{label}.ivar")
    subprocess.run(f"samtools mpileup -aa -A -d 8000 -B -Q 0 --reference {REF} {sorted_bam} "
                   f"| ivar variants -p {ivar_p} -q 20 -t 0.005 -m 10 -r {REF}",
                   shell=True, capture_output=True)

    lofreq_out = os.path.join(V3DIR, f"mix_{label}.lofreq.vcf")
    try:
        subprocess.run(f"lofreq call --no-default-filter -f {REF} -o {lofreq_out} {sorted_bam}",
                       shell=True, capture_output=True, timeout=600)
    except subprocess.TimeoutExpired:
        open(lofreq_out, "w").write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    ivar_path = f"{ivar_p}.tsv"
    ivar_all = pd.DataFrame()
    if os.path.exists(ivar_path):
        try:
            ivar_all = pd.read_csv(ivar_path, sep="\t")
            ivar_all = ivar_all[(ivar_all["REF"].str.len()==1) & (ivar_all["ALT"].str.len()==1) & (~ivar_all["ALT"].isin(["+","-"]))].copy()
            ivar_all.loc[:,"key"] = ivar_all["REGION"] + ":" + ivar_all["POS"].astype(str) + ":" + ivar_all["ALT"]
        except:
            pass

    lofreq_records = []
    if os.path.exists(lofreq_out):
        with open(lofreq_out) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8 or len(parts[3])!=1 or len(parts[4])!=1:
                    continue
                info = dict(x.split("=",1) for x in parts[7].split(";") if "=" in x)
                lofreq_records.append({"key": f"{parts[0]}:{parts[1]}:{parts[4]}", "af": float(info.get("AF",0))})

    for thresh in [0.01, 0.02, 0.03, 0.05]:
        ivar_keys = set(ivar_all[ivar_all["ALT_FREQ"] >= thresh]["key"]) if not ivar_all.empty else set()
        lofreq_keys = {r["key"] for r in lofreq_records if r["af"] >= thresh}
        concordant = ivar_keys & lofreq_keys

        tp = len(concordant & truth_keys)
        fp = len(concordant - truth_keys)
        fn = len(truth_keys - concordant)
        n_truth = len(truth_keys)
        sens = tp / n_truth if n_truth else 0
        total_neg = 13136 - n_truth
        spec = 1 - fp / total_neg if total_neg > 0 else 0

        results.append({"mixture": label, "spike_af": frac_b, "threshold": f"{thresh*100:.0f}%",
                         "tp": tp, "fn": fn, "fp": fp, "sensitivity": sens, "specificity": spec})

sep = "=" * 85
print(f"\n{sep}")
print(f"{'V3 SYNTHETIC SPIKE-IN RESULTS':^85}")
print(f"{sep}")
print(f"Truth set: {len(truth_keys)} consensus-level SNPs between samples A and B")
header = f"{'Mixture':>10} {'Spike AF':>10} {'Thresh':>8} {'TP':>5} {'FN':>5} {'FP':>5} {'Sensitivity':>12} {'Specificity':>12}"
print(header)
print("-" * 85)
for r in results:
    s_pass = "PASS" if r["sensitivity"] >= 0.80 else ""
    print(f"{r['mixture']:>10} {r['spike_af']:>9.1%} {r['threshold']:>8} {r['tp']:>5} {r['fn']:>5} {r['fp']:>5} {r['sensitivity']:>11.1%} {r['specificity']:>11.4%}  {s_pass}")

pd.DataFrame(results).to_csv(os.path.join(V3DIR, "spikein_results.tsv"), sep="\t", index=False)
print(f"\nResults saved to {V3DIR}/spikein_results.tsv")
