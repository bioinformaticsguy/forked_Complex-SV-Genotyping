#!/usr/bin/env python3
"""
Annotate a VCF with GGTyper genotype metrics from *_genotype_results.tsv(.gz).
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict


FORMAT_HEADERS = [
    '##FORMAT=<ID=GGT_GT,Number=1,Type=String,Description="GGTyper genotype call">',
    '##FORMAT=<ID=GGT_GQ,Number=1,Type=Float,Description="GGTyper mean genotype quality">',
    '##FORMAT=<ID=GGT_GQ_LOW,Number=1,Type=Float,Description="GGTyper genotype quality lower bootstrap bound">',
    '##FORMAT=<ID=GGT_GQ_HIGH,Number=1,Type=Float,Description="GGTyper genotype quality upper bootstrap bound">',
    '##FORMAT=<ID=GGT_CERT,Number=1,Type=Float,Description="GGTyper genotype certainty">',
    '##FORMAT=<ID=GGT_TR,Number=1,Type=Integer,Description="GGTyper total read pairs used">',
    '##FORMAT=<ID=GGT_OUTLIERS,Number=1,Type=Integer,Description="GGTyper outlier read pairs">',
    '##FORMAT=<ID=GGT_AVGMQ,Number=1,Type=Float,Description="GGTyper average mapping quality">',
    '##FORMAT=<ID=GGT_MINMQ,Number=1,Type=Integer,Description="GGTyper minimum mapping quality">',
    '##FORMAT=<ID=GGT_MAXMQ,Number=1,Type=Integer,Description="GGTyper maximum mapping quality">',
    '##FORMAT=<ID=GGT_QPASS,Number=1,Type=Integer,Description="GGTyper quality pass flag: 1=true, 0=false">',
    '##FORMAT=<ID=GGT_DIFF,Number=1,Type=Float,Description="GGTyper difficulty estimate">',
]

INFO_HEADERS = [
    '##INFO=<ID=GGT_VARIANT,Number=.,Type=String,Description="GGTyper variant ID(s) matched to this VCF record">',
    '##INFO=<ID=GGT_N,Number=1,Type=Integer,Description="Number of samples with GGTyper results for this record">',
    '##INFO=<ID=GGT_PASS_N,Number=1,Type=Integer,Description="Number of samples with GGTyper QualityPass=true">',
    '##INFO=<ID=GGT_MAX_GQ,Number=1,Type=Float,Description="Maximum GGTyper mean genotype quality across matched samples">',
    '##INFO=<ID=GGT_MAX_CERT,Number=1,Type=Float,Description="Maximum GGTyper certainty across matched samples">',
    '##INFO=<ID=GGT_SAMPLES,Number=.,Type=String,Description="Samples with GGTyper results for this record">',
]

FORMAT_MAP = [
    ("GGT_GT", "Genotype"),
    ("GGT_GQ", "Mean_Quality"),
    ("GGT_GQ_LOW", "Lower_Bound"),
    ("GGT_GQ_HIGH", "Upper_Bound"),
    ("GGT_CERT", "Certainty"),
    ("GGT_TR", "TotalReads"),
    ("GGT_OUTLIERS", "Outliers"),
    ("GGT_AVGMQ", "AvgMapQ"),
    ("GGT_MINMQ", "MinMapQ"),
    ("GGT_MAXMQ", "MaxMapQ"),
    ("GGT_QPASS", "QualityPass"),
    ("GGT_DIFF", "Difficulty"),
]


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def info_escape(value):
    return str(value).replace(";", "_").replace(" ", "_").replace("\t", "_").replace(",", "%2C")


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def as_quality_pass(value):
    if str(value).upper() == "TRUE":
        return "1"
    if str(value).upper() == "FALSE":
        return "0"
    return "."


def read_results(path):
    by_record_id = defaultdict(dict)
    by_variant = defaultdict(dict)
    total_rows = 0

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Variant", "Sample", "Genotype"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required GGTyper result columns: {', '.join(sorted(missing))}")

        for row in reader:
            total_rows += 1
            variant = row["Variant"]
            sample = row["Sample"]
            by_variant[variant][sample] = row
            for record_id in variant.split("|"):
                by_record_id[record_id][sample] = row

    return by_record_id, by_variant, total_rows


def add_headers(header_lines):
    existing_format = {
        line.split("ID=", 1)[1].split(",", 1)[0]
        for line in header_lines
        if line.startswith("##FORMAT=<ID=")
    }
    existing_info = {
        line.split("ID=", 1)[1].split(",", 1)[0]
        for line in header_lines
        if line.startswith("##INFO=<ID=")
    }

    additions = []
    for line in FORMAT_HEADERS:
        field_id = line.split("ID=", 1)[1].split(",", 1)[0]
        if field_id not in existing_format:
            additions.append(line)
    for line in INFO_HEADERS:
        field_id = line.split("ID=", 1)[1].split(",", 1)[0]
        if field_id not in existing_info:
            additions.append(line)

    chrom_idx = next(i for i, line in enumerate(header_lines) if line.startswith("#CHROM"))
    return header_lines[:chrom_idx] + additions + header_lines[chrom_idx:]


def append_info(info, key, value):
    entry = f"{key}={value}"
    if info in ("", "."):
        return entry
    return info + ";" + entry


def format_value(row, source_key):
    if row is None:
        return "."
    if source_key not in row or row[source_key] == "":
        return "."
    if source_key == "QualityPass":
        return as_quality_pass(row[source_key])
    return row[source_key]


def annotate_sample(format_keys, sample_value, row=None):
    values = sample_value.split(":") if sample_value not in ("", ".") else []
    data = {}
    for idx, key in enumerate(format_keys):
        data[key] = values[idx] if idx < len(values) and values[idx] != "" else "."

    for target_key, source_key in FORMAT_MAP:
        data[target_key] = format_value(row, source_key)

    new_keys = list(format_keys)
    for target_key, _ in FORMAT_MAP:
        if target_key not in new_keys:
            new_keys.append(target_key)

    return new_keys, ":".join(data.get(key, ".") for key in new_keys)


def annotate_vcf(vcf_path, results_path, output_path):
    by_record_id, by_variant, total_rows = read_results(results_path)
    matched_records = 0
    unmatched_records = 0
    annotated_sample_values = 0

    with open_text(vcf_path) as inp, open(output_path, "w") as out:
        header = []
        for line in inp:
            line = line.rstrip("\n")
            header.append(line)
            if line.startswith("#CHROM"):
                break

        header = add_headers(header)
        for line in header:
            out.write(line + "\n")

        chrom_line = header[-1].split("\t")
        sample_names = chrom_line[9:] if len(chrom_line) > 9 else []

        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                out.write(line + "\n")
                continue

            record_id = fields[2]
            sample_rows = by_record_id.get(record_id, {})

            if sample_rows:
                matched_records += 1
                variants = sorted({row["Variant"] for row in sample_rows.values()})
                samples = sorted(sample_rows)
                gqs = [as_float(row.get("Mean_Quality")) for row in sample_rows.values()]
                certs = [as_float(row.get("Certainty")) for row in sample_rows.values()]
                gqs = [value for value in gqs if value is not None]
                certs = [value for value in certs if value is not None]
                pass_n = sum(1 for row in sample_rows.values() if str(row.get("QualityPass", "")).upper() == "TRUE")

                fields[7] = append_info(fields[7], "GGT_VARIANT", ",".join(info_escape(v) for v in variants))
                fields[7] = append_info(fields[7], "GGT_N", str(len(samples)))
                fields[7] = append_info(fields[7], "GGT_PASS_N", str(pass_n))
                if gqs:
                    fields[7] = append_info(fields[7], "GGT_MAX_GQ", f"{max(gqs):.6g}")
                if certs:
                    fields[7] = append_info(fields[7], "GGT_MAX_CERT", f"{max(certs):.6g}")
                fields[7] = append_info(fields[7], "GGT_SAMPLES", ",".join(info_escape(s) for s in samples))

                if len(fields) > 8 and sample_names:
                    format_keys = fields[8].split(":") if fields[8] not in ("", ".") else []
                    updated_keys = None
                    for idx, sample in enumerate(sample_names, start=9):
                        if idx >= len(fields):
                            break
                        row = sample_rows.get(sample)
                        updated_keys, fields[idx] = annotate_sample(format_keys, fields[idx], row)
                        if row is not None:
                            annotated_sample_values += 1
                    if updated_keys is not None:
                        fields[8] = ":".join(updated_keys)
            else:
                unmatched_records += 1

            out.write("\t".join(fields) + "\n")

    print(f"GGTyper result rows read     : {total_rows}", file=sys.stderr)
    print(f"Distinct GGTyper variants    : {len(by_variant)}", file=sys.stderr)
    print(f"VCF records annotated        : {matched_records}", file=sys.stderr)
    print(f"VCF records without GGTyper  : {unmatched_records}", file=sys.stderr)
    print(f"Sample FORMAT values updated : {annotated_sample_values}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="Input VCF/VCF.gz to annotate")
    parser.add_argument("ggtyper_results", help="GGTyper *_genotype_results.tsv or .tsv.gz")
    parser.add_argument("output", help="Output annotated VCF, uncompressed")
    args = parser.parse_args()
    annotate_vcf(args.vcf, args.ggtyper_results, args.output)


if __name__ == "__main__":
    main()
