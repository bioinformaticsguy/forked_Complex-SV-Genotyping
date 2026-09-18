#!/usr/bin/env python3
"""
Generate minimal synthetic test data for the GGTyper Snakemake pipeline.

Creates two test samples in the same directory structure as real data:

  test_data/
  ├── reference/
  │   ├── test_ref.fa
  │   └── test_ref.fa.fai
  ├── TEST001_DNA_01/
  │   ├── alignment/
  │   │   ├── TEST001_DNA_01_sorted_md.bam
  │   │   └── TEST001_DNA_01_sorted_md.bam.bai
  │   └── call_sv/genome/
  │       ├── case_TEST001_DNA_01_sv.vcf.gz
  │       └── case_TEST001_DNA_01_sv.vcf.gz.tbi
  ├── TEST002_DNA_01/
      └── ...
  └── TEST003_DNA_01/
      └── ...  (151 bp reads)

Requirements:
    samtools, bcftools  (available in the genotyping conda env)

Usage:
    python3 workflow/scripts/generate_test_data.py
    # → writes test_data/ and read-length-specific sample sheets
"""

import os
import random
import subprocess
import textwrap

# ── Parameters ──────────────────────────────────────────────────────────────

OUTDIR   = "test_data"
SAMPLE_READ_LENGTHS = {
    "TEST001_DNA_01": 150,
    "TEST002_DNA_01": 150,
    "TEST003_DNA_01": 151,
}
SAMPLES  = list(SAMPLE_READ_LENGTHS)
SAMPLE_SHEETS = {
    "sample_sheets/test_samples.tsv": ["TEST001_DNA_01", "TEST002_DNA_01"],
    "sample_sheets/test_151bp_samples.tsv": ["TEST003_DNA_01"],
    "sample_sheets/test_mixed_samples.tsv": SAMPLES,
}
CHROM    = "chr1"
CHROM2   = "chr22"
CHROM_LEN = 200_000
INSERT_MEAN = 400
INSERT_SD   = 50
N_PAIRS     = 3_000      # read pairs per sample — enough for insert-size profile
SEED        = 42


# ── Helpers ─────────────────────────────────────────────────────────────────

def run(cmd, **kw):
    print(f"  $ {cmd}")
    subprocess.run(cmd, shell=True, check=True, **kw)


def make_dirs(sample):
    for d in [
        f"{OUTDIR}/{sample}/alignment",
        f"{OUTDIR}/{sample}/call_sv/genome",
    ]:
        os.makedirs(d, exist_ok=True)


# ── Reference ───────────────────────────────────────────────────────────────

def make_reference():
    ref_dir = f"{OUTDIR}/reference"
    os.makedirs(ref_dir, exist_ok=True)
    ref_path = f"{ref_dir}/test_ref.fa"

    random.seed(SEED)
    bases = "ACGT"

    with open(ref_path, "w") as f:
        for chrom, length in [(CHROM, CHROM_LEN), (CHROM2, CHROM_LEN)]:
            seq = "".join(random.choices(bases, k=length))
            f.write(f">{chrom}\n")
            for i in range(0, length, 80):
                f.write(seq[i:i+80] + "\n")

    run(f"samtools faidx {ref_path}")
    print(f"  Reference: {ref_path}")
    return ref_path


# ── BAM ─────────────────────────────────────────────────────────────────────

def make_bam(sample, read_len, seed_offset=0):
    random.seed(SEED + seed_offset)

    sam_path = f"{OUTDIR}/{sample}/alignment/temp.sam"
    bam_base = f"{OUTDIR}/{sample}/alignment/{sample}_sorted_md"

    reads = []   # list of (pos1, pos2, sam_line_R1, sam_line_R2)

    for i in range(N_PAIRS):
        insert = max(250, int(random.gauss(INSERT_MEAN, INSERT_SD)))
        pos1   = random.randint(1, CHROM_LEN - insert - read_len)
        pos2   = pos1 + insert - read_len

        name = f"read_{i:06d}"
        seq  = "".join(random.choices("ACGT", k=read_len))
        qual = "I" * read_len

        r1 = (f"{name}\t99\t{CHROM}\t{pos1}\t60\t{read_len}M"
              f"\t=\t{pos2}\t{insert}\t{seq}\t{qual}"
              f"\tRG:Z:RG1\n")
        r2 = (f"{name}\t147\t{CHROM}\t{pos2}\t60\t{read_len}M"
              f"\t=\t{pos1}\t{-insert}\t{seq}\t{qual}"
              f"\tRG:Z:RG1\n")

        reads.append((pos1, pos2, r1, r2))

    # Write coordinate-sorted SAM (sort by first-in-pair position)
    reads.sort(key=lambda x: x[0])

    with open(sam_path, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:{CHROM}\tLN:{CHROM_LEN}\n")
        f.write(f"@SQ\tSN:{CHROM2}\tLN:{CHROM_LEN}\n")
        f.write(f"@RG\tID:RG1\tSM:{sample}\tPL:ILLUMINA\tLB:lib1\n")
        for _, _, r1, r2 in reads:
            f.write(r1)
            f.write(r2)

    # Re-sort properly with samtools (handles R2 ordering)
    run(f"samtools sort -o {bam_base}.bam {sam_path}")
    run(f"samtools index {bam_base}.bam")
    os.remove(sam_path)

    print(f"  BAM: {bam_base}.bam")
    return f"{bam_base}.bam"


# ── VCF ─────────────────────────────────────────────────────────────────────

def svdb_info(variant_id, chrom, pos, end, svtype, svlen, caller,
              imprecise=False, mateid=None):
    """Build a minimal svdb-merged INFO field with only declared header tags."""
    parts = [f"END={end}", f"SVTYPE={svtype}", f"SVLEN={svlen}"]
    if imprecise:
        parts.append("IMPRECISE")
    if mateid:
        parts.append(f"MATEID={mateid}")
    parts += [f"set={caller}", "FOUNDBY=1", f"svdb_origin={caller}", "SUPP_VEC=010"]
    return ";".join(parts)


def make_vcf(sample, seed_offset=0):
    vcf_out = f"{OUTDIR}/{sample}/call_sv/genome/case_{sample}_sv.vcf.gz"

    # Fixed SVs on standard chromosomes (same for both samples so they merge)
    variants = [
        # (ID, CHROM, POS, REF, ALT, QUAL, END, SVTYPE, SVLEN, CALLER, IMPRECISE, GT)
        # ── Manta PRECISE calls ──────────────────────────────────────────────
        ("MantaDEL:1:0:0:0:0:0",       CHROM,  20000, "A", "<DEL>",       200,  25000, "DEL",       -5000,  "manta", False, "0/1"),
        ("MantaDUP:TANDEM:1:0:0:0:0:0",CHROM,  50000, "C", "<DUP:TANDEM>",180,  52000, "DUP",        2000,  "manta", False, "0/1"),
        ("MantaINV:1:0:0:0:0:0",       CHROM,  70000, "T", "<INV>",       150,  80000, "INV",       10000,  "manta", False, "0/1"),
        # ── Manta BND pair ───────────────────────────────────────────────────
        ("MantaBND:1:0:0:0:0:0:0",     CHROM,  90000, "N", f"N[{CHROM}:120000[", 120, 90000, "BND", 0, "manta", False, "0/1"),
        ("MantaBND:1:0:0:0:0:0:1",     CHROM, 120000, "N", f"]{CHROM}:90000]N",  120,120000, "BND", 0, "manta", False, "0/1"),
        # ── Manta IMPRECISE (should be filtered out by converter) ────────────
        ("MantaDEL:2:0:0:0:0:0",       CHROM, 140000, "G", "<DEL>",        80, 160000, "DEL",      -20000, "manta", True,  "0/1"),
        # ── CNVnator (should be filtered by --caller manta) ──────────────────
        ("CNVnator_del_1",              CHROM,  30000, "N", "<DEL>",         0,  35000, "DEL",       -5000, "cnvnator", False, "1/1"),
        # ── TIDDIT (filtered by --caller manta, but present in merged VCF) ──
        ("SV_1_1",                      CHROM,  45000, "N", "<DEL>",       100,  47000, "DEL",       -2000, "tiddit",   False, "0/1"),
        # ── Small Manta INS (filtered by converter as unsupported) ───────────
        ("MantaINS:1:0:0:0:0:0",        CHROM,  60000, "A",
         "AGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", 300, 60000, "INS", 50, "manta", False, "1/1"),
    ]

    vcf_lines = []
    vcf_lines.append("##fileformat=VCFv4.1")
    vcf_lines.append(f"##FILTER=<ID=PASS,Description=\"All filters passed\">")
    vcf_lines.append(f"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">")
    vcf_lines.append(f"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">")
    vcf_lines.append(f"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">")
    vcf_lines.append(f"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise SV\">")
    vcf_lines.append(f"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"BND mate ID\">")
    vcf_lines.append(f"##INFO=<ID=svdb_origin,Number=1,Type=String,Description=\"Caller origin\">")
    vcf_lines.append(f"##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Support vector\">")
    vcf_lines.append(f"##INFO=<ID=set,Number=1,Type=String,Description=\"Caller set\">")
    vcf_lines.append(f"##INFO=<ID=FOUNDBY,Number=1,Type=Integer,Description=\"Found by N callers\">")
    vcf_lines.append(f"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    vcf_lines.append(f"##contig=<ID={CHROM},length={CHROM_LEN}>")
    vcf_lines.append(f"##contig=<ID={CHROM2},length={CHROM_LEN}>")
    vcf_lines.append(f"##svdbcmdline=svdb --merge --pass_only --vcf ...")
    vcf_lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")

    for (vid, chrom, pos, ref, alt, qual, end,
         svtype, svlen, caller, imprecise, gt) in variants:

        mateid = None
        if svtype == "BND" and vid.endswith(":0"):
            mateid = vid[:-1] + "1"
        elif svtype == "BND" and vid.endswith(":1"):
            mateid = vid[:-1] + "0"

        info = svdb_info(vid, chrom, pos, end, svtype, svlen,
                         caller, imprecise, mateid)
        vcf_lines.append(
            f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}"
            f"\tPASS\t{info}\tGT\t{gt}"
        )

    vcf_text = "\n".join(vcf_lines) + "\n"

    tmp = vcf_out.replace(".vcf.gz", "_tmp.vcf")
    with open(tmp, "w") as f:
        f.write(vcf_text)

    run(f"bcftools sort -O z -o {vcf_out} {tmp}")
    run(f"bcftools index -t {vcf_out}")
    os.remove(tmp)

    print(f"  VCF: {vcf_out}")
    return vcf_out


# ── Sample sheet ─────────────────────────────────────────────────────────────

def make_sample_sheet(tsv_path, samples):
    os.makedirs("sample_sheets", exist_ok=True)

    with open(tsv_path, "w") as f:
        f.write("sample_id\tbam_path\tvcf_path\n")
        for sample in samples:
            bam = os.path.abspath(f"{OUTDIR}/{sample}/alignment/{sample}_sorted_md.bam")
            vcf = os.path.abspath(f"{OUTDIR}/{sample}/call_sv/genome/case_{sample}_sv.vcf.gz")
            f.write(f"{sample}\t{bam}\t{vcf}\n")

    print(f"  Sample sheet: {tsv_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=== Generating test data ===")

    print(f"\n[1/{len(SAMPLES)+2}] Reference genome")
    make_reference()

    for i, sample in enumerate(SAMPLES):
        print(f"\n[{i+2}/{len(SAMPLES)+2}] Sample: {sample}")
        make_dirs(sample)
        make_bam(sample, SAMPLE_READ_LENGTHS[sample], seed_offset=i * 1000)
        make_vcf(sample, seed_offset=i * 1000)

    print(f"\n[{len(SAMPLES)+2}/{len(SAMPLES)+2}] Sample sheets")
    for tsv_path, samples in SAMPLE_SHEETS.items():
        make_sample_sheet(tsv_path, samples)

    print("\n=== Done ===")
    print(f"Test data written to: {os.path.abspath(OUTDIR)}/")
    print(f"Run the pipeline with:")
    print(f"  snakemake --snakefile workflow/Snakefile --configfile configs/config.yaml \\")
    print(f"    --config samples_sheet=sample_sheets/test_samples.tsv output_dir=output -n")


if __name__ == "__main__":
    main()
