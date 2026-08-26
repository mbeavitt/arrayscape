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
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument("--period", type=int,
                     help="take every class of this period, e.g. 178 -> 178_*")
    sel.add_argument("--class", dest="klass", help="one exact class name")
    ap.add_argument("--genome", required=True, help="name to put in the headers")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--max-len", type=int, default=0, help="0 = no limit")
    args = ap.parse_args()

    csv.field_size_limit(1 << 30)
    # TRASH's class suffix is a per-run index, not a stable identity: the same
    # CEN178 array is 178_1 in one assembly and 178_2 or 178_436 in the next,
    # so selecting one exact name silently loses whole genomes. Match on the
    # period instead, which is the part that means something.
    prefix = f"{args.period}_" if args.period else None
    n = 0
    used = set()
    with open(args.csv, newline="") as fh:
        for r in csv.DictReader(fh):
            cls = r.get("class") or ""
            if prefix is not None:
                if not cls.startswith(prefix):
                    continue
            elif cls != args.klass:
                continue
            used.add(cls)
            seq = (r.get("sequence") or "").strip()
            if not seq or len(seq) < args.min_len:
                continue
            if args.max_len and len(seq) > args.max_len:
                continue
            print(f">{args.genome}|{r['seqID']}:{r['start']}-{r['end']}")
            print(seq.upper())
            n += 1
    want = f"period {args.period}" if args.period else f"class {args.klass}"
    print(f"{args.csv}: {n:,} monomers of {want}"
          + (f" [{', '.join(sorted(used))}]" if used else ""), file=sys.stderr)
    if not n:
        # a genome with none of this repeat is a fact about the genome, not an
        # error: warn and leave an empty file so one assembly cannot sink a run
        print(f"WARNING: {args.csv}: no rows of {want}", file=sys.stderr)


if __name__ == "__main__":
    main()
