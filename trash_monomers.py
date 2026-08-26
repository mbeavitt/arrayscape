#!/usr/bin/env python3
"""TRASH repeats table -> a monomer fasta for one genome.

trash-py writes <name>_repeats_with_seq.csv with a row per repeat unit and its
sequence, already in scaffold coordinates, so unlike the FasTAN route there is
nothing to rebase -- this only selects a class and relabels the headers into
the <genome>|<chrom>:<start>-<end> form the rest of the pipeline reads.

    trash_monomers.py --class 178_1 --genome Col-0 table.csv > Col-0.fasta
"""

import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="*_repeats_with_seq.csv")
    ap.add_argument("--class", dest="klass", required=True, help="repeat class")
    ap.add_argument("--genome", required=True, help="name to put in the headers")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--max-len", type=int, default=0, help="0 = no limit")
    args = ap.parse_args()

    csv.field_size_limit(1 << 30)
    n = 0
    with open(args.csv, newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("class") != args.klass:
                continue
            seq = (r.get("sequence") or "").strip()
            if not seq or len(seq) < args.min_len:
                continue
            if args.max_len and len(seq) > args.max_len:
                continue
            print(f">{args.genome}|{r['seqID']}:{r['start']}-{r['end']}")
            print(seq.upper())
            n += 1
    print(f"{args.csv}: {n:,} monomers of class {args.klass}", file=sys.stderr)
    if not n:
        sys.exit(f"{args.csv}: no rows of class {args.klass}")


if __name__ == "__main__":
    main()
