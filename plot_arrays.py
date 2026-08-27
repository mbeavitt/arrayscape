#!/usr/bin/env python3
"""Draw monomer arrays from a chromomer table: one row per genome, one plot
per chromosome.

Reads the table `chromomer` writes and pulls the genome, chromosome and
position out of the id column, so a single table covering many assemblies
becomes one figure per chromosome with the assemblies stacked as rows.

    chromomer monomers.fasta > colours.tsv
    plot_arrays.py colours.tsv -o plots/

Ids are parsed with --id-regex, which must supply named groups `genome`,
`chrom` and `start` (`end` optional). The default matches

    <genome>|<chrom>:<start>-<end>
"""

import argparse
import os
import re
import sys

import numpy as np

from chromomer import colours, scaling

DEFAULT_ID_RE = r"(?P<genome>[^|]+)\|(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)"
EMPTY = 0.955                      # a bin of chromosome holding no monomer
MIN_SHARED = 10                    # bins two genomes must share to be compared


def read_table(path, id_re):
    """id column -> genome, chromosome, midpoint; plus the pc columns."""
    pat = re.compile(id_re)
    genomes, chroms, mids, pcs = [], [], [], []
    skipped = 0
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        need = ("id", "pc1", "pc2", "pc3")
        missing = [c for c in need if c not in col]
        if missing:
            sys.exit(f"{path}: missing column(s) {', '.join(missing)}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = pat.match(f[col["id"]])
            if not m:
                skipped += 1
                continue
            d = m.groupdict()
            genomes.append(d["genome"])
            chroms.append(d["chrom"])
            start = int(d["start"])
            end = int(d.get("end") or start)
            mids.append((start + end) // 2)
            pcs.append((float(f[col["pc1"]]), float(f[col["pc2"]]),
                        float(f[col["pc3"]])))
    if not genomes:
        sys.exit("no ids matched --id-regex; check the id format")
    if skipped:
        # silence here would be dangerous: a regex that fits only part of the
        # table would quietly plot a subset and look perfectly fine
        print(f"warning: {skipped:,} of {skipped + len(genomes):,} ids did not "
              f"match --id-regex and were dropped", file=sys.stderr)
    gnames, gidx = np.unique(genomes, return_inverse=True)
    cnames, cidx = np.unique(chroms, return_inverse=True)
    return (gnames, gidx, cnames, cidx, np.array(mids, np.int64),
            np.array(pcs, np.float32))


def tracks(gidx, mids, pcs, ngenomes, nbins, align, smooth, keep=0.99):
    """Mean position in composition space per (genome, bin), and the window.

    Averaging the embedding and colouring afterwards, rather than averaging
    colours: mixing colours pulls every bin towards grey, whereas mixing
    coordinates keeps a bin of similar monomers as saturated as one monomer.

    With --align each genome is shifted by the median position of its own
    monomers, so centromeres line up even when the arms differ in length.
    """
    pos = mids.astype(np.float64)
    if align:
        shift = np.zeros(ngenomes)
        for g in np.unique(gidx):
            shift[g] = np.median(pos[gidx == g])
        pos = pos - shift[gidx]
    tail = (1 - keep) / 2 * 100
    lo, hi = np.percentile(pos, [tail, 100 - tail])
    pad = 0.15 * (hi - lo)                    # context either side of the arrays
    lo, hi = lo - pad, hi + pad
    if hi <= lo:                              # every monomer at one position
        lo, hi = lo - 1, hi + 1

    b = np.clip(((pos - lo) / (hi - lo) * nbins).astype(int), 0, nbins - 1)
    flat = gidx.astype(np.int64) * nbins + b
    n = np.bincount(flat, minlength=ngenomes * nbins).reshape(ngenomes, nbins)
    sums = np.stack([np.bincount(flat, weights=pcs[:, c],
                                 minlength=ngenomes * nbins).reshape(ngenomes, nbins)
                     for c in range(3)], axis=2)

    occupied = n > 0
    if smooth > 1:
        # Rolling mean over neighbouring columns: at this zoom a column holds
        # only ~10 monomers, so raw column means are mostly sampling noise.
        kern = np.ones(smooth)
        roll = lambda a: np.apply_along_axis(
            lambda v: np.convolve(v, kern, "same"), 1, a)
        sums = np.stack([roll(sums[:, :, c]) for c in range(3)], axis=2)
        denom = roll(n.astype(float))
    else:
        denom = n.astype(float)
    mean = np.where(occupied[:, :, None],
                    sums / np.maximum(denom, 1)[:, :, None], 0).astype(np.float32)
    return mean, n, lo, hi


def render(mean, n):
    """Colour the filled bins, scaling over the bins actually drawn.

    Each chromosome is normalised against its own median bin, so colour is
    departure from this chromosome's consensus and is not comparable between
    figures.
    """
    img = np.full(mean.shape, EMPTY, np.float32)
    filled = n > 0
    if filled.any():
        img[filled] = colours(mean[filled].astype(np.float64))
    return img


def distances(mean, n):
    """Mean distance in composition space between each pair of genomes.

    Compared only over the bins both genomes fill, so a genome with a short
    array is judged on the array it has rather than on the emptiness around it.
    """
    g = len(mean)
    filled = n > 0
    d = np.full((g, g), np.nan)
    for i in range(g):
        for j in range(i + 1, g):
            both = filled[i] & filled[j]
            if both.sum() < MIN_SHARED:
                continue
            v = np.linalg.norm(mean[i, both] - mean[j, both], axis=1).mean()
            d[i, j] = d[j, i] = v
    return d


def order_from_file(path, names):
    """Row order taken from a file, one genome per line.

    Anything the file does not mention keeps its relative order at the end,
    so a partial list still works and nothing is silently dropped.
    """
    want = [l.strip() for l in open(path) if l.strip()]
    rank = {g: i for i, g in enumerate(want)}
    return np.array(sorted(range(len(names)),
                           key=lambda i: (rank.get(names[i], len(rank)), names[i])))


def order_rows(mean, n, how, names):
    """Row order: alphabetical, or each genome next to the one it matches best.

    Walks a nearest-neighbour chain: start at the most atypical genome, then
    repeatedly step to whichever genome not yet placed is closest to the one
    just placed. Neighbouring rows are therefore the most similar pairs, which
    a one-dimensional projection does not guarantee.
    """
    if how == "name":
        return np.argsort(names)
    d = distances(mean, n)
    g = len(mean)
    # a genome sharing too few bins with everything has an all-NaN row, so
    # average only over the rows that have any comparison at all
    comparable = ~np.all(np.isnan(d), axis=1)
    typical = np.full(g, -np.inf)
    if comparable.any():
        typical[comparable] = np.nanmean(d[comparable], axis=1)
    start = int(np.argmax(typical))

    order, left = [start], set(range(g)) - {start}
    while left:
        row = d[order[-1]]
        cand = [j for j in left if not np.isnan(row[j])]
        # nothing comparable left: restart the chain at whatever remains
        nxt = min(cand, key=lambda j: row[j]) if cand else min(left)
        order.append(nxt)
        left.discard(nxt)
    return np.array(order)


def colour_key(ax, mean, n, res=160):
    """The colour field itself over the PC1/PC2 plane, on the plot's own scale.

    A reference, not a data display: it shows what any given colour means,
    which a scatter of the points cannot, since the points hide the gamut
    wherever they happen not to fall.
    """
    v = mean[n > 0].astype(np.float64)
    centre, scale = scaling(v)
    gx, gy = np.meshgrid(np.linspace(-1, 1, res), np.linspace(-1, 1, res))
    grid = centre + np.stack([gx.ravel() * scale[0], gy.ravel() * scale[1],
                              np.zeros(gx.size)], axis=1)
    ax.imshow(colours(grid, ref=v).reshape(res, res, 3), origin="lower",
              extent=(-1, 1, -1, 1), interpolation="bilinear")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("PC1", fontsize=7, labelpad=1)
    ax.set_ylabel("PC2", fontsize=7, labelpad=1)
    for s in ax.spines.values():
        s.set_color("#c9c8c3")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="tsv written by chromomer")
    ap.add_argument("-o", "--outdir", default=".", help="where to write the pngs")
    ap.add_argument("--align", action="store_true",
                    help="centre each genome on its own median monomer position")
    ap.add_argument("--bins", type=int, default=3000, help="horizontal resolution")
    ap.add_argument("--smooth", type=int, default=9, metavar="N",
                    help="rolling mean over N columns (default 9, 1 = off)")
    ap.add_argument("--order", choices=("name", "similarity"), default="similarity",
                    help="row order (default similarity)")
    ap.add_argument("--order-from", metavar="FILE",
                    help="take the row order from this file, one genome per line "
                         "(overrides --order; the order used is always written "
                         "alongside each figure as <chrom>.order.txt)")
    ap.add_argument("--id-regex", default=DEFAULT_ID_RE,
                    help="regex with named groups genome, chrom, start[, end]")
    ap.add_argument("--no-key", action="store_true", help="omit the PC1/PC2 key")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gnames, gidx, cnames, cidx, mids, pcs = read_table(args.table, args.id_regex)
    os.makedirs(args.outdir, exist_ok=True)
    print(f"{len(mids):,} monomers  {len(gnames)} genomes  "
          f"{len(cnames)} chromosomes", file=sys.stderr)

    for ci, chrom in enumerate(cnames):
        sel = cidx == ci
        mean, n, lo, hi = tracks(gidx[sel], mids[sel], pcs[sel], len(gnames),
                                 args.bins, args.align, args.smooth)
        keep = n.sum(1) > 0                       # genomes present on this chromosome
        mean, n = mean[keep], n[keep]
        names = gnames[keep]
        img = render(mean, n)
        row = (order_from_file(args.order_from, names) if args.order_from
               else order_rows(mean, n, args.order, names))
        img, mean, n, names = img[row], mean[row], n[row], names[row]

        fig, ax = plt.subplots(figsize=(14, 11), constrained_layout=True)
        ax.imshow(img, interpolation="nearest", aspect="auto",
                  extent=(lo / 1e6, hi / 1e6, len(names), 0))
        ax.set_yticks(np.arange(len(names)) + 0.5)
        ax.set_yticklabels(names, fontsize=3.1)
        ax.tick_params(axis="y", length=0, pad=1)
        ax.set_xlabel("offset from centromere (Mb)" if args.align else "position (Mb)",
                      fontsize=10, color="#52514e")
        ax.set_title(chrom, fontsize=12, loc="left", pad=8, color="#0b0b0b")
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        if not args.no_key:
            k = fig.add_axes((0.885, 0.80, 0.075, 0.075 * 14 / 11))
            colour_key(k, mean, n)
        with open(os.path.join(args.outdir, f"{chrom}.order.txt"), "w") as fh:
            fh.write("\n".join(names) + "\n")
        out = os.path.join(args.outdir, f"{chrom}.png")
        fig.savefig(out, dpi=260, facecolor="white")
        plt.close(fig)
        print(f"  {out}  {len(names)} genomes", file=sys.stderr)


if __name__ == "__main__":
    main()
