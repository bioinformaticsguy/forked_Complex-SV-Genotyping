#!/usr/bin/env python3
"""
Filter a GGTyper-annotated VCF using FORMAT/GGT_* genotype evidence.
"""
import argparse
import gzip
import re
import sys


MISSING = {"", "."}


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_float(value):
    if value in MISSING:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value):
    if value in MISSING:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def genotype_has_variant(gt):
    if gt in MISSING:
        return False
    alleles = [allele for allele in re.split(r"[\/|]", gt) if allele]
    if not alleles:
        return False
    return any(allele not in {"REF", "."} for allele in alleles)


def sample_data(format_keys, sample_value):
    values = sample_value.split(":") if sample_value not in MISSING else []
    return {
        key: values[idx] if idx < len(values) and values[idx] != "" else "."
        for idx, key in enumerate(format_keys)
    }


def sample_passes(data, args):
    if not genotype_has_variant(data.get("GGT_GT", ".")):
        return False

    if args.require_quality_pass and data.get("GGT_QPASS", ".") != "1":
        return False

    checks = [
        (parse_float(data.get("GGT_CERT", ".")), args.min_certainty),
        (parse_float(data.get("GGT_GQ", ".")), args.min_genotype_quality),
        (parse_int(data.get("GGT_TR", ".")), args.min_total_reads),
        (parse_float(data.get("GGT_AVGMQ", ".")), args.min_avg_mapq),
    ]
    for observed, threshold in checks:
        if threshold is not None and (observed is None or observed < threshold):
            return False

    return True


def add_filter_metadata(header_lines, args):
    criteria = [
        f"require_quality_pass={int(args.require_quality_pass)}",
        f"min_certainty={args.min_certainty}",
        f"min_genotype_quality={args.min_genotype_quality}",
        f"min_total_reads={args.min_total_reads}",
        f"min_avg_mapq={args.min_avg_mapq}",
        f"min_passing_samples={args.min_passing_samples}",
    ]
    line = '##GGTyperFilter=<Description="Post-genotyping GGTyper FORMAT filter",Criteria="{}">'.format(
        ",".join(criteria)
    )
    chrom_idx = next(i for i, line in enumerate(header_lines) if line.startswith("#CHROM"))
    return header_lines[:chrom_idx] + [line] + header_lines[chrom_idx:]


def filter_vcf(input_vcf, output_vcf, args):
    total = kept = missing_format = 0
    passing_samples_total = 0

    with open_text(input_vcf) as inp, open(output_vcf, "w") as out:
        header = []
        for line in inp:
            line = line.rstrip("\n")
            header.append(line)
            if line.startswith("#CHROM"):
                break

        header = add_filter_metadata(header, args)
        for line in header:
            out.write(line + "\n")

        chrom_fields = header[-1].split("\t")
        sample_names = chrom_fields[9:] if len(chrom_fields) > 9 else []

        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue

            total += 1
            fields = line.split("\t")
            if len(fields) < 10 or not sample_names:
                missing_format += 1
                continue

            format_keys = fields[8].split(":") if fields[8] not in MISSING else []
            required_keys = {"GGT_GT", "GGT_QPASS", "GGT_CERT", "GGT_GQ", "GGT_TR", "GGT_AVGMQ"}
            if not required_keys.issubset(format_keys):
                missing_format += 1
                continue

            passing = 0
            for sample_idx in range(9, min(len(fields), 9 + len(sample_names))):
                if sample_passes(sample_data(format_keys, fields[sample_idx]), args):
                    passing += 1

            passing_samples_total += passing
            if passing >= args.min_passing_samples:
                kept += 1
                out.write(line + "\n")

    print(f"VCF records read              : {total}", file=sys.stderr)
    print(f"VCF records kept              : {kept}", file=sys.stderr)
    print(f"VCF records removed           : {total - kept}", file=sys.stderr)
    print(f"Records missing GGT FORMAT    : {missing_format}", file=sys.stderr)
    print(f"Passing sample calls observed : {passing_samples_total}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_vcf", help="Input GGTyper-annotated VCF/VCF.gz")
    parser.add_argument("output_vcf", help="Output filtered VCF, uncompressed")
    parser.add_argument("--min-certainty", type=float, default=0.8)
    parser.add_argument("--min-genotype-quality", type=float, default=20)
    parser.add_argument("--min-total-reads", type=int, default=10)
    parser.add_argument("--min-avg-mapq", type=float, default=20)
    parser.add_argument("--min-passing-samples", type=int, default=1)
    quality_pass = parser.add_mutually_exclusive_group()
    quality_pass.add_argument(
        "--require-quality-pass",
        dest="require_quality_pass",
        action="store_true",
        help="Require FORMAT/GGT_QPASS=1 for a sample call to pass",
    )
    quality_pass.add_argument(
        "--no-require-quality-pass",
        dest="require_quality_pass",
        action="store_false",
        help="Do not require FORMAT/GGT_QPASS=1 for a sample call to pass",
    )
    parser.set_defaults(require_quality_pass=True)
    args = parser.parse_args()

    filter_vcf(args.input_vcf, args.output_vcf, args)


if __name__ == "__main__":
    main()
