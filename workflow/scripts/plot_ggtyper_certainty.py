#!/usr/bin/env python3
"""
Plot the GGTyper Certainty score distribution from *_genotype_results.tsv(.gz).
"""
import argparse
import csv
import gzip
import math
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_bool(value):
    return str(value).strip().upper() == "TRUE"


def has_variant(genotype):
    if genotype in ("", "."):
        return False
    alleles = genotype.replace("|", "/").split("/")
    return any(allele not in ("REF", ".") for allele in alleles)


def read_certainties(path):
    rows = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Certainty", "Genotype", "QualityPass"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            value = row.get("Certainty", "")
            try:
                certainty = float(value)
            except ValueError:
                continue
            if not math.isfinite(certainty):
                continue
            rows.append(
                {
                    "certainty": certainty,
                    "genotype": row.get("Genotype", ""),
                    "quality_pass": parse_bool(row.get("QualityPass", "")),
                    "variant_call": has_variant(row.get("Genotype", "")),
                }
            )
    return rows


def subset_rows(rows, mode):
    if mode == "all":
        return rows
    if mode == "variant":
        return [row for row in rows if row["variant_call"]]
    if mode == "quality-pass":
        return [row for row in rows if row["quality_pass"]]
    if mode == "variant-quality-pass":
        return [row for row in rows if row["variant_call"] and row["quality_pass"]]
    raise ValueError(mode)


def quantile(values, q):
    if not values:
        return None
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - pos) + ordered[hi] * (pos - lo)


def histogram(values, bins):
    counts = [0] * bins
    for value in values:
        idx = int(value * bins)
        if idx >= bins:
            idx = bins - 1
        if idx < 0:
            idx = 0
        counts[idx] += 1
    return counts


def write_summary(path, rows, thresholds):
    groups = [
        ("all", rows),
        ("variant", subset_rows(rows, "variant")),
        ("quality-pass", subset_rows(rows, "quality-pass")),
        ("variant-quality-pass", subset_rows(rows, "variant-quality-pass")),
    ]
    quantiles = [0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0]

    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "n", "mean"] + [f"q{int(q * 100):02d}" for q in quantiles])
        for name, group_rows in groups:
            values = [row["certainty"] for row in group_rows]
            mean = sum(values) / len(values) if values else None
            writer.writerow(
                [name, len(values), format_float(mean)]
                + [format_float(quantile(values, q)) for q in quantiles]
            )

    threshold_path = path.with_suffix(".thresholds.tsv")
    with open(threshold_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group", "threshold", "kept", "removed", "kept_fraction"])
        for name, group_rows in groups:
            values = [row["certainty"] for row in group_rows]
            for threshold in thresholds:
                kept = sum(1 for value in values if value >= threshold)
                removed = len(values) - kept
                frac = kept / len(values) if values else 0
                writer.writerow([name, threshold, kept, removed, f"{frac:.6f}"])


def format_float(value):
    if value is None:
        return "NA"
    return f"{value:.6f}"


def svg_text(x, y, text, size=12, anchor="start", weight="normal", color="#1f2937"):
    return (
        f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" '
        f'font-family="Arial, sans-serif" font-weight="{weight}" fill="{color}">{text}</text>'
    )


def write_svg(path, values, thresholds, title, bins=50):
    width, height = 980, 580
    left, right, top, bottom = 82, 32, 72, 82
    plot_w = width - left - right
    plot_h = height - top - bottom
    counts = histogram(values, bins)
    max_count = max(counts) if counts else 1
    nice_max = max(1, math.ceil(max_count / 10) * 10)
    bar_gap = 2
    bar_w = plot_w / bins

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(left, 34, title, 22, weight="bold"),
        svg_text(left, 56, f"n = {len(values):,}", 13, color="#4b5563"),
    ]

    for i in range(6):
        y_value = nice_max * i / 5
        y = top + plot_h - (y_value / nice_max) * plot_h
        parts.append(f'<line x1="{left}" x2="{width - right}" y1="{y:.2f}" y2="{y:.2f}" stroke="#e5e7eb"/>')
        parts.append(svg_text(left - 10, y + 4, str(int(y_value)), 11, anchor="end", color="#6b7280"))

    for i, count in enumerate(counts):
        x = left + i * bar_w
        h = (count / nice_max) * plot_h
        y = top + plot_h - h
        parts.append(
            f'<rect x="{x + bar_gap / 2:.2f}" y="{y:.2f}" width="{max(0.5, bar_w - bar_gap):.2f}" '
            f'height="{h:.2f}" fill="#2563eb"/>'
        )

    for threshold in thresholds:
        x = left + threshold * plot_w
        parts.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top}" y2="{top + plot_h}" stroke="#dc2626" stroke-width="2" stroke-dasharray="6 5"/>')
        parts.append(svg_text(x + 5, top + 18, str(threshold), 12, color="#991b1b"))

    parts.extend(
        [
            f'<line x1="{left}" x2="{width - right}" y1="{top + plot_h}" y2="{top + plot_h}" stroke="#111827"/>',
            f'<line x1="{left}" x2="{left}" y1="{top}" y2="{top + plot_h}" stroke="#111827"/>',
            svg_text(width / 2, height - 28, "GGTyper certainty score", 14, anchor="middle"),
            svg_text(22, top + plot_h / 2, "Number of calls", 14, anchor="middle"),
        ]
    )

    for tick in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        x = left + tick * plot_w
        parts.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{top + plot_h}" y2="{top + plot_h + 6}" stroke="#111827"/>')
        parts.append(svg_text(x, top + plot_h + 24, f"{tick:.1f}", 11, anchor="middle", color="#4b5563"))

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="GGTyper *_genotype_results.tsv or .tsv.gz")
    parser.add_argument("--output-prefix", default=None, help="Output path prefix")
    parser.add_argument(
        "--mode",
        choices=["all", "variant", "quality-pass", "variant-quality-pass"],
        default="all",
        help="Calls to include in the histogram",
    )
    parser.add_argument("--bins", type=int, default=50)
    parser.add_argument(
        "--thresholds",
        default="0.5,0.7,0.8,0.9,0.95",
        help="Comma-separated threshold lines and retention summary values",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    prefix = Path(args.output_prefix) if args.output_prefix else input_path.with_suffix("").with_suffix("")
    thresholds = [float(value) for value in args.thresholds.split(",") if value.strip()]
    rows = read_certainties(input_path)
    selected = subset_rows(rows, args.mode)
    values = [row["certainty"] for row in selected]

    if not values:
        raise SystemExit(f"No rows selected for mode '{args.mode}'")

    svg_path = prefix.parent / f"{prefix.name}.certainty_{args.mode}.svg"
    summary_path = prefix.parent / f"{prefix.name}.certainty_summary.tsv"
    write_svg(svg_path, values, thresholds, f"GGTyper certainty distribution ({args.mode})", args.bins)
    write_summary(summary_path, rows, thresholds)

    print(f"Wrote plot: {svg_path}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote threshold table: {summary_path.with_suffix('.thresholds.tsv')}")


if __name__ == "__main__":
    main()
