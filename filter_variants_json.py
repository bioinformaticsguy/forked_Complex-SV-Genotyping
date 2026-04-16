#!/usr/bin/env python3
"""
Filter variants.json to only keep variants where all breakpoints
are on standard chromosomes (chr1-22, chrX, chrY, chrM).
"""
import json
import argparse

STANDARD = {f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y', 'M']}


def is_standard(variant):
    for allele in variant.get('VAR', {}).values():
        for junction in allele.values():
            if junction.get('rNameLeft',  '') not in STANDARD:
                return False
            if junction.get('rNameRight', '') not in STANDARD:
                return False
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',  help='Input variants.json')
    p.add_argument('output', help='Output filtered variants.json')
    args = p.parse_args()

    variants = json.load(open(args.input))
    filtered = {k: v for k, v in variants.items() if is_standard(v)}

    json.dump(filtered, open(args.output, 'w'), indent=4)
    print(f'Kept    : {len(filtered)} variants')
    print(f'Removed : {len(variants) - len(filtered)} variants on non-standard contigs')


if __name__ == '__main__':
    main()
