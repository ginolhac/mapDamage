#!/usr/bin/env python

import mapdamage
import itertools
import csv
import subprocess
import sys

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
    path_to_seqtk = mapdamage.rscript.construct_path("seqtk",folder="seqtk") 
    bases = {"A":0,"C":0,"G":0,"T":0}
    alp = ["A","C","G","T"]
    try:
        proc = subprocess.Popen([path_to_seqtk,"comp",filename],stdout=subprocess.PIPE)
        out = proc.communicate()[0]
        for li in out.splitlines():
            tabs = li.split() # 1 is the ref, 2 is the total and then the base counts A, C, G and T.
            bases["A"] = bases["A"] + int(tabs[2])
            bases["C"] = bases["C"] + int(tabs[3])
            bases["G"] = bases["G"] + int(tabs[4])
            bases["T"] = bases["T"] + int(tabs[5])
    except OSError or ValueError:
        sys.stderr.write("Error: Seqtk failed\n")
        sys.exit(1)
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
