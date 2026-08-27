#!/usr/bin/env python3
"""Put every monomer in one rotational frame, and one orientation.

An annotator picks where a repeat unit starts, and that choice is arbitrary and
made per array: FasTAN chooses a partition phase per array from its own seed, so
two arrays of the same satellite come out cut at different points in the repeat.
The sequences are the same, but their k-mer composition is not, and the shift is
constant across every monomer of an array, so it does not average out -- it
becomes a whole-array offset in composition space.

Measured on 155 A. thaliana assemblies, that offset dominated: 73% of the
variance along the leading component sat *between* centromeres rather than
within them, against 16% for an annotator that phases consistently. Comparing
arrays or genomes without rephasing largely compares phase choices.

Each array is rotated onto a common reference monomer, forwards or reverse
complemented, whichever matches better -- so this fixes orientation at the same
time, and two annotators become comparable to each other as well.

    rephase.py --build-reference monomers.fasta > ref.fasta
    rephase.py --reference ref.fasta monomers.fasta > phased.fasta
"""

import argparse
import re
import sys

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
ID_RE = re.compile(r"(?P<genome>[^|]+)\|(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)")
CONS_SAMPLE = 40                  # monomers per array used to build its consensus


def read_fasta(path):
    ids, seqs, name, parts = [], [], None, []
    fh = sys.stdin.buffer if path == "-" else open(path, "rb")
    with fh:
        for line in fh:
            if line[:1] == b">":
                if name is not None:
                    ids.append(name)
                    seqs.append(b"".join(parts))
                name, parts = line[1:].split()[0].decode(), []
            elif name is not None:
                parts.append(line.strip())
    if name is not None:
        ids.append(name)
        seqs.append(b"".join(parts))
    return ids, seqs


def arrays(ids, unit, gap):
    """Group monomers into arrays: same genome and chromosome, and no positional
    jump bigger than `gap`. Phase is constant within an array, not between."""
    keyed = []
    for i, name in enumerate(ids):
        m = ID_RE.match(name)
        keyed.append((m["genome"], m["chrom"], int(m["start"]), i) if m
                     else (name, "", 0, i))
    keyed.sort()
    out, cur, prev = [], [], None
    for g, c, s, i in keyed:
        if prev and (g, c) == (prev[0], prev[1]) and s - prev[2] <= gap:
            cur.append(i)
        else:
            if cur:
                out.append(cur)
            cur = [i]
        prev = (g, c, s)
    if cur:
        out.append(cur)
    return out


def consensus(monomers, unit):
    """Column-wise consensus of the unit-length monomers, else the first one."""
    full = [m for m in monomers if len(m) == unit][:CONS_SAMPLE]
    if not full:
        m = next((m for m in monomers if len(m) >= unit), None)
        return np.frombuffer(m[:unit], np.uint8) if m else None
    stack = np.vstack([np.frombuffer(m, np.uint8) for m in full])
    counts = np.stack([(stack == b).sum(axis=0) for b in BASES])
    return BASES[counts.argmax(axis=0)]


def rotations(cons, ref, unit):
    """Best (offset, score) over every rotation of cons against ref."""
    doubled = np.concatenate([cons, cons])
    scores = np.array([(doubled[r:r + unit] == ref).sum() for r in range(unit)])
    r = int(scores.argmax())
    return r, int(scores[r])


def best_rotation(cons, ref, unit):
    """(shift, identity, strand) for an array.

    Arrays occur in both orientations, so the reference is matched against the
    consensus and its reverse complement; the better of the two wins and brings
    the orientation with it.
    """
    fwd_r, fwd_s = rotations(cons, ref, unit)
    rc = np.frombuffer(cons.tobytes().translate(COMPLEMENT)[::-1], np.uint8)
    rev_r, rev_s = rotations(rc, ref, unit)
    if rev_s > fwd_s:
        return rev_r, rev_s, "-"
    return fwd_r, fwd_s, "+"


def apply(seq, rot, strand):
    """Rotate one monomer into the array's frame, flipping it first if needed."""
    if strand == "-":
        seq = seq.translate(COMPLEMENT)[::-1]
    if not seq:
        return seq
    r = rot % len(seq)
    return seq[r:] + seq[:r]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="monomer fasta; `-` reads stdin")
    ap.add_argument("--reference", help="reference monomer fasta")
    ap.add_argument("--build-reference", action="store_true",
                    help="emit a reference monomer instead: the consensus of the "
                         "largest array, an arbitrary but fixed phase")
    ap.add_argument("--period", type=int, default=178, help="repeat unit (default 178)")
    ap.add_argument("--gap", type=int, default=0, metavar="BP",
                    help="positional jump that starts a new array "
                         "(default 10x --period)")
    ap.add_argument("--min-identity", type=float, default=0.0,
                    help="skip arrays matching the reference below this fraction")
    args = ap.parse_args()
    unit = args.period
    gap = args.gap or unit * 10

    ids, seqs = read_fasta(args.fasta)
    if not ids:
        sys.exit("no sequences found")
    groups = arrays(ids, unit, gap)

    if args.build_reference:
        big = max(groups, key=len)
        cons = consensus([seqs[i] for i in big], unit)
        if cons is None:
            sys.exit("largest array has no usable monomer")
        m = ID_RE.match(ids[big[0]])
        where = f"{m['genome']} {m['chrom']}:{m['start']}" if m else ids[big[0]]
        print(f">rephase_ref {where} n={len(big)}")
        print(cons.tobytes().decode())
        return

    if not args.reference:
        sys.exit("need --reference REF.fasta (or --build-reference)")
    _, rseq = read_fasta(args.reference)
    ref = np.frombuffer(rseq[0][:unit], np.uint8)
    if len(ref) != unit:
        sys.exit(f"reference is {len(ref)} bp, expected {unit}")

    flipped = skipped = 0
    out = sys.stdout
    for grp in groups:
        cons = consensus([seqs[i] for i in grp], unit)
        if cons is None:
            skipped += len(grp)
            continue
        rot, ident, strand = best_rotation(cons, ref, unit)
        if ident / unit < args.min_identity:
            skipped += len(grp)
            continue
        flipped += len(grp) if strand == "-" else 0
        for i in grp:
            out.write(f">{ids[i]}\n{apply(seqs[i], rot, strand).decode()}\n")

    print(f"{len(groups):,} arrays, {len(ids):,} monomers; "
          f"{flipped:,} reverse complemented, {skipped:,} skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
