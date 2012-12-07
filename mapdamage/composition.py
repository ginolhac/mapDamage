#!/usr/bin/env python

import mapdamage
import itertools
import csv


def count_ref_comp(read, chrom, before, after, comp):
  """ record basae composition in external genomic regions """
  std = '-' if read.is_reverse else '+'

  _update_table(comp[chrom]['5p'][std], before, xrange(-len(before), 0))
  _update_table(comp[chrom]['3p'][std], after,  xrange(1, len(after) + 1))


def count_read_comp(read, chrom, length, comp):
  """ record base composition of read, discard marked nucleotides """
  std, seq = '+', read.query
  if read.is_reverse:
    std, seq = '-', mapdamage.seq.revcomp(seq)

  _update_table(comp[chrom]['5p'][std], seq,           xrange(1, length + 1))
  _update_table(comp[chrom]['3p'][std], reversed(seq), xrange(-1, - length - 1, -1))


def _update_table(table, sequence, indices):
  for (index, nt) in itertools.izip(indices, sequence):
    if nt in table:
      table[nt][index] += 1


def get_base_comp(filename,destination=False):
    """
    Gets the basecomposition of all the sequences in filename
    and returns the value to destination if given.
    """
    f = open(filename,"r")
    bases = {"A":0,"C":0,"G":0,"T":0}
    for li in f:
        if li[0] == ">":
            continue
        for b in li:
            b = b.upper()
            if (b in bases.keys()):
                bases[b] = bases[b] + 1 
    f.close()
    su = sum(bases.values())
    for i in bases.keys():
        bases[i] = float(bases[i])/float(su)
    if (destination==False):
        return bases
    else:
        # write the results
        fo = open(destination,"w")
        fow = csv.DictWriter(fo,bases.keys())
        fow.writeheader()
        fow.writerow(bases)
        fo.close
