import mapdamage
import itertools
import csv
import subprocess
import sys

import mapdamage.seqtk as seqtk


def count_ref_comp(read, chrom, before, after, comp):
  """ record basae composition in external genomic regions """
  std = '-' if read.is_reverse else '+'

  _update_table(comp[chrom]['5p'][std], before, range(-len(before), 0))
  _update_table(comp[chrom]['3p'][std], after,  range(1, len(after) + 1))


def count_read_comp(read, chrom, length, comp):
  """ record base composition of read, discard marked nucleotides """
  std, seq = '+', read.query
  if read.is_reverse:
    std, seq = '-', mapdamage.seq.revcomp(seq)

  _update_table(comp[chrom]['5p'][std], seq,           range(1, length + 1))
  _update_table(comp[chrom]['3p'][std], reversed(seq), range(-1, - length - 1, -1))


def _update_table(table, sequence, indices):
  for (index, nt) in zip(indices, sequence):
    if nt in table:
      table[nt][index] += 1


def get_base_comp(filename,destination=False):
    """
    Gets the basecomposition of all the sequences in filename
    and returns the value to destination if given.
    """
    bases = {"A":0,"C":0,"G":0,"T":0}
    alp = ["A","C","G","T"]

    for stats in seqtk.comp(filename).values():
        for key in bases:
            bases[key] += stats[key]

    # get the base frequencies
    ba_su = sum(bases.values())
    for ba in alp:
        bases[ba] = float(bases[ba])/float(ba_su)
    if (destination==False):
        return bases
    else:
        # write the results
        fo = open(destination,"w")
        vals = [str(bases[i]) for i in alp]
        fo.write(",".join(alp)+"\n")
        fo.write(",".join(vals)+"\n")
        fo.close()

def read_base_comp(filename):
    """
    Read the base compition from a file created by get_base_comp
    """
    fh = open(filename)
    first_line = True
    for li in fh:
        li = li.rstrip()
        lp = li.split()
        if first_line:
            header = lp
            first_line = False
        else:
            body = lp
    bases = {}
    for ba,per in zip(header,body):
        bases[ba] = per
    return bases
