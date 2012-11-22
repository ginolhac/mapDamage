#!/usr/bin/env python

import mapdamage
import itertools


def count_ref_comp(read, chrom, before, after, comp):
  std = '-' if read.is_reverse else '+'

  _update_table(comp[chrom]['5p'][std], before, xrange(-len(before), 0))
  _update_table(comp[chrom]['3p'][std], after,  xrange(1, len(after) + 1))


def count_read_comp(read, chrom, length, comp):
  std, seq = '+', read.query
  if read.is_reverse:
    std, seq = '-', mapdamage.seq.revcomp(seq)

  _update_table(comp[chrom]['5p'][std], seq,           xrange(1, length + 1))
  _update_table(comp[chrom]['3p'][std], reversed(seq), xrange(-1, - length - 1, -1))


def count_read_comp_with_qual(read, chrom, opt, comp):
  std, seq, qual = '+', read.query, read.qqual
  if read.is_reverse:
    std, seq, qual = '-', mapdamage.seq.revcomp(seq), qual[::-1]

  _update_table_with_qual(comp[chrom]['5p'][std], seq, qual, \
      xrange(1, opt.length + 1), opt.minqual)
  _update_table_with_qual(comp[chrom]['3p'][std], reversed(seq), reversed(qual), \
      xrange(-1, - opt.length - 1, -1), opt.minqual)


def _update_table(table, sequence, indices):
  for (index, nt) in itertools.izip(indices, sequence):
    if nt in table:
      table[nt][index] += 1


def _update_table_with_qual(table, sequence, qualities, indices, threshold):
  for (index, nt, qual) in itertools.izip(indices, sequence, qualities):
    if nt in table and (ord(qual)-33) >= threshold:
      table[nt][index] += 1
