#!/usr/bin/env python3
"""FasTAN annotations -> a monomer bed in scaffold coordinates.

Why this exists
---------------
FAtoGDB breaks every scaffold at each run of non-ACGT, and FasTAN's masks are
stored per contig, so they are in *contig* coordinates. Upstream ANOtoBED
prints those numbers beside the scaffold name, which is short by the summed
length of the preceding contigs and gaps for anything past the first gap in a
scaffold -- a median shift of about 4 Mb on a gappy A. thaliana assembly, and
one that no bounds check catches, since a shifted coordinate is still inside
the scaffold.

Rather than raise FAtoGDB's minimum gap length to stop it splitting (which
changes the annotation, not just the coordinates), this reads the contig layout
straight out of the GDB with `GDBshow -h`

    >Chr1 <0,5538260] :: Contig 1 <0,5538260>
    >Chr1 [5538360,29508191> :: Contig 2 <0,23969831>

and adds each contig's scaffold start back onto its intervals, which is what a
patched ANOtoBED does at the print. A patched ANOtoBED is detected and left
alone.

    ano_monomers.py --gdb g.1gdb --bed g.178.bed --period 178 > g.monomers.bed
"""

import argparse
import re
import subprocess
import sys

CONTIG_RE = re.compile(
    r"^>(?P<scaf>\S+)\s+[<\[](?P<sbeg>\d+),(?P<send>\d+)[\]>]\s*::\s*"
    r"Contig\s+(?P<idx>\d+)\s+<0,(?P<clen>\d+)>")


def contig_table(gdb, gdbshow="GDBshow"):
    """[(scaffold, sbeg, clen)] in contig order, from the GDB itself."""
    try:
        out = subprocess.run([gdbshow, "-h", gdb], capture_output=True,
                             text=True, check=True).stdout
    except FileNotFoundError:
        sys.exit(f"{gdbshow} not on PATH; needed to translate contig coordinates")
    except subprocess.CalledProcessError as e:
        sys.exit(f"{gdbshow} failed on {gdb}: {e.stderr.strip()}")
    rows = []
    for line in out.splitlines():
        m = CONTIG_RE.match(line)
        if m:
            rows.append((m["scaf"], int(m["sbeg"]), int(m["clen"])))
    if not rows:
        sys.exit(f"could not read a contig table from {gdbshow} -h {gdb}")
    return rows


def read_bed(path):
    """ANOtoBED output: comment lines carry the parse points for the interval
    that precedes them, so the two travel together.

    Data lines are split on TAB, not on whitespace: ANOtoBED prints the whole
    FASTA header as the first field, and outside toy assemblies that header has
    spaces in it ("CM026974.1 Rattus norvegicus strain ... chromosome 1"), so
    whitespace splitting silently reads a word as the start coordinate. Only the
    first token is kept as the sequence name, which is what samtools and every
    other tool downstream uses.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                f = line.split()
                if out and len(f) > 1 and f[1] == "Parse:":
                    out[-1][3] = [int(x) for x in f[2:]]
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            out.append([f[0].split()[0], int(f[1]), int(f[2]), None, f[3:]])
    return out


def already_scaffold_coords(bed, contigs):
    """True if the bed is already rebased (a patched ANOtoBED).

    Contig coordinates restart at zero at every contig boundary, so within a
    scaffold they are not monotonic; scaffold coordinates are. Only decide on
    a scaffold that actually has more than one contig, otherwise the two are
    the same thing and the question does not arise.
    """
    multi = {s for s, _, _ in contigs
             if sum(1 for t, _, _ in contigs if t == s) > 1}
    seen, monotonic = {}, True
    any_multi = False
    for scaf, beg, _, _, _ in bed:
        if scaf not in multi:
            continue
        any_multi = True
        if scaf in seen and beg < seen[scaf]:
            monotonic = False
        seen[scaf] = beg
    return monotonic if any_multi else True


def translate(bed, contigs):
    """Add each contig's scaffold start onto its own intervals.

    Intervals arrive grouped by contig in contig order, so a step to the next
    contig of a scaffold is marked by the coordinate going backwards, or by an
    interval that will not fit inside the contig we think we are in.
    """
    by_scaf = {}
    for scaf, sbeg, clen in contigs:
        by_scaf.setdefault(scaf, []).append((sbeg, clen))

    at, prev = {}, {}
    out, dropped = [], 0
    for scaf, beg, end, parse, rest in bed:
        ctgs = by_scaf.get(scaf)
        if not ctgs:
            dropped += 1
            continue
        i = at.get(scaf, 0)
        if scaf in prev and beg < prev[scaf]:
            i += 1
        while i < len(ctgs) - 1 and end > ctgs[i][1]:
            i += 1
        at[scaf], prev[scaf] = i, beg
        sbeg = ctgs[i][0]
        # only the interval moves: parse points are offsets from its start,
        # which is how rephase.py reads them (beg + points[i]), so shifting
        # them too would double-count the contig offset
        out.append([scaf, beg + sbeg, end + sbeg, parse, rest])
    if dropped:
        print(f"warning: {dropped} interval(s) on scaffolds absent from the GDB",
              file=sys.stderr)
    return out


def monomers(bed, period, tol):
    """Consecutive parse points inside an interval become monomers.

    Parse points are offsets from the interval start, not coordinates, so a
    monomer runs from beg+points[i] to beg+points[i+1].
    """
    lo, hi = period * (1 - tol), period * (1 + tol)
    for scaf, beg, end, parse, rest in bed:
        if not parse or len(parse) < 2:
            continue                      # need two boundaries to make a monomer
        strand = "-" if any(r == "-" for r in rest) else "+"
        for a, b in zip(parse, parse[1:]):
            if lo <= b - a <= hi:
                yield scaf, beg + a, beg + b, strand


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gdb", help="the .1gdb the .1ano came from; only needed "
                                  "when contig coordinates must be translated")
    ap.add_argument("--bed", required=True, help="ANOtoBED output for one class")
    ap.add_argument("--period", type=int, required=True, help="repeat period, e.g. 178")
    ap.add_argument("--tol", type=float, default=0.15,
                    help="fractional length tolerance on a monomer (default 0.15)")
    ap.add_argument("--gdbshow", default="GDBshow")
    ap.add_argument("--assume", choices=("auto", "contig", "scaffold"), default="auto",
                    help="coordinate space of the bed (default auto-detect)")
    args = ap.parse_args()

    bed = read_bed(args.bed)
    if not bed:
        sys.exit(f"{args.bed}: no intervals")

    space = args.assume
    if space == "auto":
        if not args.gdb:
            sys.exit("--gdb is needed to detect the coordinate space; "
                     "pass --assume scaffold or --assume contig if you know it")
        space = ("scaffold" if already_scaffold_coords(bed, contig_table(args.gdb, args.gdbshow))
                 else "contig")
    contigs = contig_table(args.gdb, args.gdbshow) if space == "contig" else []
    if space == "contig":
        bed = translate(bed, contigs)
        print(f"{args.bed}: translated contig -> scaffold coordinates "
              f"({len(contigs)} contigs)", file=sys.stderr)
    else:
        print(f"{args.bed}: already in scaffold coordinates, left alone",
              file=sys.stderr)

    n = 0
    for scaf, a, b, strand in monomers(bed, args.period, args.tol):
        print(f"{scaf}\t{a}\t{b}\t.\t0\t{strand}")
        n += 1
    print(f"{args.bed}: {n:,} monomers", file=sys.stderr)


if __name__ == "__main__":
    main()
