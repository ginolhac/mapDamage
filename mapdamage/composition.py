#!/usr/bin/env python

import mapdamage
import itertools


def countRefComp(read, chrom, before, after, comp):
  std = '-' if read.is_reverse else '+'

  _update_table(comp[chrom]['5p'][std], before, xrange(-len(before), 0))
  _update_table(comp[chrom]['3p'][std], after,  xrange(1, len(after) + 1))


def countReadComp(read, chrom, lg, comp):
  std, seq = '+', read.query
  if read.is_reverse:
    std, seq = '-', mapdamage.seq.revcomp(seq)

  _update_table(comp[chrom]['5p'][std], seq,           xrange(1, lg + 1))
  _update_table(comp[chrom]['3p'][std], reversed(seq), xrange(-1, - lg - 1, -1))


def _update_table(table, sequence, indices):
  for (index, nt) in itertools.izip(indices, sequence):
    if nt in table:
      table[nt][index]      += 1
