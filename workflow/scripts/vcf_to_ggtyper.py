#!/usr/bin/env python3
"""
Convert SV VCF (Manta calls) to GGTyper JSON format.
Handles DEL, DUP, INV, and BND records. Skips IMPRECISE calls.

BND ALT notation → GGTyper directions:
  N[chr:pos[   dirLeft=left,  dirRight=right  (deletion-like)
  N]chr:pos]   dirLeft=left,  dirRight=left   (inversion left junction)
  [chr:pos[N   dirLeft=right, dirRight=right  (inversion right junction)
  ]chr:pos]N   dirLeft=right, dirRight=left   (dup-like)

Usage:
    python vcf_to_ggtyper.py input.vcf.gz output.json [--caller manta] [--max N]
"""
import re
import json
import gzip
import hashlib
import argparse

_MAX_VID = 200  # profile filename = vid + ".profile"; Linux limit is 255 bytes

def safe_vid(vid):
    """Return vid unchanged if short enough, otherwise a stable SHA-256 prefix."""
    if len(vid) <= _MAX_VID:
        return vid
    return "sv_" + hashlib.sha256(vid.encode()).hexdigest()[:20]


def parse_info(info_str):
    info = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            info[k] = v
        else:
            info[field] = True
    return info


def parse_bnd_alt(alt):
    """
    Parse BND ALT field.
    Returns (mate_chrom, mate_pos, dir_left, dir_right) or None.
    dir_left/right = direction at current POS / mate POS respectively.
    """
    # N[chr:pos[  → left of current, right of mate
    m = re.match(r'^[A-Za-z.]+\[([^:[\]]+):(\d+)\[$', alt)
    if m:
        return m.group(1), int(m.group(2)), "left", "right"
    # N]chr:pos]  → left of current, left of mate
    m = re.match(r'^[A-Za-z.]+\]([^:[\]]+):(\d+)\]$', alt)
    if m:
        return m.group(1), int(m.group(2)), "left", "left"
    # [chr:pos[N  → right of current, right of mate
    m = re.match(r'^\[([^:[\]]+):(\d+)\[[A-Za-z.]+$', alt)
    if m:
        return m.group(1), int(m.group(2)), "right", "right"
    # ]chr:pos]N  → right of current, left of mate
    m = re.match(r'^\]([^:[\]]+):(\d+)\][A-Za-z.]+$', alt)
    if m:
        return m.group(1), int(m.group(2)), "right", "left"
    return None


def make_junction(chrom_l, pos_l, dir_l, chrom_r, pos_r, dir_r):
    return {
        "rNameLeft":      chrom_l,
        "xLeft":          pos_l,
        "directionLeft":  dir_l,
        "rNameRight":     chrom_r,
        "xRight":         pos_r,
        "directionRight": dir_r,
    }


def convert(vcf_path, output_path, caller_filter=None, skip_imprecise=True, max_variants=None):
    variants = {}
    skipped_imprecise = 0
    skipped_caller    = 0
    skipped_other     = 0
    skipped_bnd_mate  = 0
    count = 0

    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos   = int(fields[1])
            vid   = fields[2]
            alt   = fields[4]
            info  = parse_info(fields[7])

            # Filter by caller
            if caller_filter:
                origin = info.get('svdb_origin', '')
                if caller_filter not in origin:
                    skipped_caller += 1
                    continue

            # Skip imprecise
            if skip_imprecise and 'IMPRECISE' in info:
                skipped_imprecise += 1
                continue

            svtype = info.get('SVTYPE', '')
            end    = int(info.get('END', pos))

            junctions = None

            if svtype == 'DEL':
                junctions = {
                    "1": make_junction(chrom, pos, "left", chrom, end, "right")
                }

            elif svtype == 'DUP':
                # Tandem dup: LEFT of END joins RIGHT of POS
                junctions = {
                    "1": make_junction(chrom, end, "left", chrom, pos, "right")
                }

            elif svtype == 'INV':
                junctions = {
                    "1": make_junction(chrom, pos, "left",  chrom, end, "left"),
                    "2": make_junction(chrom, pos, "right", chrom, end, "right"),
                }

            elif svtype == 'BND':
                # Only process the first mate of each pair to avoid duplicates.
                # Manta names mates ...X:0 and ...X:1; skip the :1 mate.
                mateid = info.get('MATEID', '')
                if vid > mateid:
                    skipped_bnd_mate += 1
                    continue

                parsed = parse_bnd_alt(alt)
                if parsed is None:
                    skipped_other += 1
                    continue

                mate_chrom, mate_pos, dir_l, dir_r = parsed
                junctions = {
                    "1": make_junction(chrom, pos, dir_l, mate_chrom, mate_pos, dir_r)
                }

            else:
                skipped_other += 1
                continue

            # Group junctions under the affected chromosome(s)
            chrom_key = chrom
            variants[safe_vid(vid)] = {
                "VAR": {
                    chrom_key: junctions
                }
            }
            count += 1

            if max_variants and count >= max_variants:
                break

    with open(output_path, 'w') as f:
        json.dump(variants, f, indent=4)

    print(f"Converted  : {count} variants -> {output_path}")
    print(f"Skipped    : {skipped_imprecise} IMPRECISE, "
          f"{skipped_caller} wrong caller, "
          f"{skipped_bnd_mate} BND second-mates, "
          f"{skipped_other} unsupported (INS/etc.)")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('vcf',    help='Input VCF or VCF.gz')
    p.add_argument('output', help='Output JSON path')
    p.add_argument('--caller', default='manta',
                   help='Filter to svdb_origin value (default: manta)')
    p.add_argument('--all-callers', action='store_true',
                   help='Do not filter by caller')
    p.add_argument('--keep-imprecise', action='store_true',
                   help='Include IMPRECISE calls (not recommended)')
    p.add_argument('--max', type=int, default=None,
                   help='Maximum number of variants to convert')
    args = p.parse_args()

    caller = None if args.all_callers else args.caller
    convert(
        args.vcf,
        args.output,
        caller_filter=caller,
        skip_imprecise=not args.keep_imprecise,
        max_variants=args.max,
    )


if __name__ == '__main__':
    main()
