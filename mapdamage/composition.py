#!/usr/bin/env python

import mapdamage


def countRefComp(read, chrom, before, after, comp):
  std = '-' if read.is_reverse else '+'

  for (i, nt) in enumerate(before.upper()):
    if nt in mapdamage.seq.letters:
      # record base composition in the reference
      j = i - len(before)
      comp[chrom]['5p'][std]['Total'][j] += 1
      comp[chrom]['5p'][std][nt][j] += 1

  for (i, nt) in enumerate(after.upper()):
    if nt in mapdamage.seq.letters:
      comp[chrom]['3p'][std]['Total'][i+1] += 1
      comp[chrom]['3p'][std][nt][i+1] += 1


def countReadComp(read, chrom, lg, comp):
  seq = read.query
  std = '+'
  if read.is_reverse:
    std = '-'
    seq = mapdamage.seq.revcomp(seq)

  # compute base composition for base aligned only
  for (i, nt) in enumerate(seq.upper()):
    if nt in mapdamage.seq.letters and i < lg:
      comp[chrom]['5p'][std]['Total'][i+1] += 1
      comp[chrom]['5p'][std][nt][i+1] += 1

  rev = seq[::-1]
  for (i, nt) in enumerate(rev.upper()):
    if nt in mapdamage.seq.letters and i < lg:
      j = (i * -1)
      j -= 1
      # record base composition
      comp[chrom]['3p'][std]['Total'][j] += 1
      comp[chrom]['3p'][std][nt][j] += 1



