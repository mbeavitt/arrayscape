# arrayscape

Draw tandem-repeat arrays across many genomes: annotate every assembly, cut
the monomers out, colour them with [chromomer](https://github.com/mbeavitt/chromomer),
and draw one figure per chromosome with the assemblies stacked as rows.

The embedding is fitted once over every monomer from every genome, which is the
whole point — a colour then means the same thing in every row and every
chromosome, so a block of shared colour is a real statement about shared
composition rather than an artefact of each assembly being scaled on its own.

## Install

Needs `chromomer`, GNU `parallel`, and whichever annotator you use:

```
pip install git+https://github.com/mbeavitt/chromomer
```

- **TRASH** route (default): `trash-py` on PATH.
- **FasTAN** route: `FAtoGDB`, `FasTAN`, `ANOtoBED`, `GDBshow` and `samtools`.

## Use

```
arrayscape -o out *.fasta                       # TRASH
arrayscape -o out --annotator fastan *.fasta    # FasTAN
```

Or start part-way in, if the slow stage is already done:

```
arrayscape -o out --from-annotations *_repeats_with_seq.csv   # after trash-py
arrayscape -o out --annotator fastan --from-annotations *.1ano
arrayscape -o out --from-monomers *.monomers.fasta
```

Output:

```
out/
  monomers/<genome>.fasta   per-genome monomers
  monomers.fasta            everything, concatenated
  chromosomes.txt           the identifier report
  colours.tsv               the chromomer table
  plots/<chrom>.png         one figure per chromosome
  logs/                     per-genome annotator and coordinate logs
```

## Chromosome identifiers

Before anything is plotted the run reports how many unique chromosome
identifiers it found and how many genomes each appears in, and warns about:

- identifiers present in fewer than half the genomes — usually unplaced
  scaffolds, or one assembly naming things differently;
- identifiers differing only in case or punctuation (`Chr1` vs `chr1`);
- headers that did not parse at all.

None of these stop the run. They matter because each variant becomes its own
band of rows in the figure, which is easy to miss and easy to misread.

## The FasTAN coordinate problem

FAtoGDB breaks every scaffold at each run of non-ACGT, and FasTAN stores its
masks per contig, so they are in **contig** coordinates. Upstream `ANOtoBED`
prints those numbers beside the *scaffold* name, so any interval past the first
gap in a scaffold is short by the summed length of the preceding contigs and
gaps — a shift of megabases on a gappy assembly, and one that no bounds check
catches, because a shifted coordinate is still inside the scaffold.

`ano_monomers.py` reads the real contig layout out of the GDB with `GDBshow -h`

```
>Chr1 <0,5538260] :: Contig 1 <0,5538260>
>Chr1 [5538360,29508191> :: Contig 2 <0,23969831>
```

and adds each contig's scaffold start back onto its own intervals — the same
correction a patched ANOtoBED makes at the print. Contig coordinates restart at
zero at every boundary while scaffold coordinates do not, so which of the two
you have is detectable, and a patched ANOtoBED is left alone. Override the
guess with `--assume contig|scaffold` if you ever need to.

This deliberately does **not** raise FAtoGDB's minimum gap length to stop it
splitting: that changes which annotation you get, not merely how it is
reported.

Note that only the interval moves. The `# Parse:` points ANOtoBED prints are
offsets from the interval's own start, not coordinates -- which is how
`rephase.py` reads them (`beg + points[i]`) -- so a monomer runs from
`beg + points[i]` to `beg + points[i+1]` and the offsets themselves must be
left alone. Shifting them as well double-counts the contig offset.

Verified by round trip on *A. thaliana* Chr1 (2 contigs, 4,738 intervals):
subtracting each contig's `sbeg` from patched output reproduces what upstream
emits, and translating that back recovers the patched intervals exactly and
yields all 15,011 monomers identically, correcting shifts of up to 5,538,360 bp.
On a whole assembly the resulting Chr1 monomers span 14,280,540-18,437,782 --
the centromere -- and come to 69,842 against 69,894 from the established
recipe, a 0.07% difference from the monomer length filter alone.

## Options

```
  -o DIR              output directory (default out)
  -j N                parallel jobs (default: nproc)
  --annotator NAME    trash (default) or fastan
  --class NAME        TRASH repeat class (default 178_1)
  --period N          repeat period for FasTAN, and the monomer length (default 178)
  --tol F             fractional length tolerance on a monomer (default 0.15)
  --from-annotations  inputs are existing TRASH csv or FasTAN .1ano
  --from-monomers     inputs are already monomer fastas named <genome>.*
  --bins N            horizontal resolution of a row (default 3000)
  --smooth N          rolling mean over N columns (default 9)
  --no-align          plot absolute position instead of offset from centromere
  --keep-temp         keep the per-genome working files
```

## Licence

MIT.
