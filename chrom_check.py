#!/usr/bin/env python3
"""Report the chromosome identifiers in a monomer fasta and flag odd ones.

A figure across a pangenome only makes sense if the assemblies agree on what to call
a chromosome. They frequently do not: one assembly says `Chr1`, another `chr1`,
a third carries unplaced scaffolds through. Each variant becomes its own row
band in the figure and the mismatch is easy to miss, so it is worth saying out
loud before anything is plotted.

    chrom_check.py monomers.fasta [--min-genomes-frac 0.5]
"""

import argparse
import collections
import re
import sys

DEFAULT_ID_RE = r"(?P<genome>[^|]+)\|(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--id-regex", default=DEFAULT_ID_RE)
    ap.add_argument("--min-genomes-frac", type=float, default=0.5,
                    help="warn below this share of genomes (default 0.5)")
    args = ap.parse_args()

    pat = re.compile(args.id_regex)
    genomes_of = collections.defaultdict(set)
    counts = collections.Counter()
    genomes = set()
    unparsed = 0
    with open(args.fasta) as fh:
        for line in fh:
            if line[0] != ">":
                continue
            m = pat.match(line[1:].split()[0])
            if not m:
                unparsed += 1
                continue
            d = m.groupdict()
            genomes.add(d["genome"])
            genomes_of[d["chrom"]].add(d["genome"])
            counts[d["chrom"]] += 1

    ng = len(genomes)
    print(f"{ng} genomes, {len(counts)} unique chromosome identifiers")
    for chrom in sorted(counts, key=lambda c: (-len(genomes_of[c]), c)):
        print(f"  {chrom:<24} {len(genomes_of[chrom]):>4}/{ng} genomes  "
              f"{counts[chrom]:>10,} monomers")

    cut = max(1, int(args.min_genomes_frac * ng))
    odd = [c for c in counts if len(genomes_of[c]) < cut]
    problems = 0
    if unparsed:
        print(f"WARNING: {unparsed:,} header(s) did not match the id pattern",
              file=sys.stderr)
        problems += 1
    if odd:
        print(f"WARNING: {len(odd)} identifier(s) present in fewer than "
              f"{cut}/{ng} genomes -- likely unplaced scaffolds or a naming "
              f"mismatch, and each becomes its own band of rows:", file=sys.stderr)
        for c in sorted(odd, key=lambda c: -counts[c])[:20]:
            print(f"    {c}  ({len(genomes_of[c])} genome(s), "
                  f"{counts[c]:,} monomers)", file=sys.stderr)
        if len(odd) > 20:
            print(f"    ... and {len(odd) - 20} more", file=sys.stderr)
        problems += 1

    # case/prefix variants of the same name are the usual cause, so name them
    norm = collections.defaultdict(list)
    for c in counts:
        norm[re.sub(r"[^0-9a-z]", "", c.lower())].append(c)
    clashes = {k: v for k, v in norm.items() if len(v) > 1}
    if clashes:
        print("WARNING: identifiers differing only in case or punctuation:",
              file=sys.stderr)
        for v in clashes.values():
            print(f"    {' / '.join(sorted(v))}", file=sys.stderr)
        problems += 1

    return 0 if problems == 0 else 0     # advisory only; never blocks the run


if __name__ == "__main__":
    sys.exit(main())
