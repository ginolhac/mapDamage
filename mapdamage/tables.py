#!/usr/bin/env python

import os
import sys
import collections

import mapdamage
from mapdamage.version import __version__


class recursivedefaultdict(collections.defaultdict):
  def __init__(self):
    self.default_factory = type(self) 


def initializeMut(ref, lg):  
  tab = recursivedefaultdict() 
  for contig in ref:
    for end in ('5p','3p'):
      for std in ('+','-'):
        for mut in mapdamage.seq.header:
          tab[contig][end][std][mut] = collections.defaultdict(int)
  
  return tab


def printMut(mut, op, out):
  out.write("# table produced by mapDamage version %s\n" % __version__)
  out.write("# using mapped file %s and %s as reference file\n" % (op.filename, op.ref))
  out.write("# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads\n")
  out.write("Chr\tEnd\tStd\tPos\t%s\t%s\n" % ("\t".join(mapdamage.seq.letters),"\t".join(mapdamage.seq.mutations)))
  for ref in sorted(mut):
    for end in mut[ref]:
      for std in mut[ref][end]:
        for i in range(op.length):
          out.write("%s\t%s\t%s\t%d" % (ref, end, std, i+1))
          for mis in mapdamage.seq.header:
            out.write("\t%d" % mut[ref][end][std][mis][i])
          out.write("\n")


def initializeComp(ref, around,lg): 
  tab = recursivedefaultdict() 
  for contig in ref:
    for end in ('5p','3p'):
      for std in ('+','-'):
        for letters in mapdamage.seq.letters:
          tab[contig][end][std][letters] = collections.defaultdict(int)

  return tab


def printComp(comp, op, out):
  out.write("# table produced by mapDamage version %s\n" % __version__)
  out.write("# using mapped file %s and %s as reference file\n" % (op.filename, op.ref))
  out.write("# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads\n")
  out.write("Chr\tEnd\tStd\tPos\t%s\n" % ("\t".join(mapdamage.seq.letters)))
  for ref in sorted(comp):
    for end in comp[ref]:
      start_i, end_i = op.around, op.length
      if end == '3p':
        start_i, end_i = end_i, start_i
        
      for std in comp[ref][end]:
        for current_i in range(-1 * start_i, end_i + 1):
          if current_i:
            out.write("%s\t%s\t%s\t%d" % (ref, end, std, current_i))
            for base in mapdamage.seq.letters:
              out.write("\t%d" % comp[ref][end][std][base][current_i])
            out.write("\n")


def initializeLg():
  tab = recursivedefaultdict() 
  for std in ('+','-'):
    tab[std] = collections.defaultdict(int)

  return tab 


def printLg(tab, op, out):
  out.write("# table produced by mapDamage version %s\n" % __version__)
  out.write("# using mapped file %s and %s as reference file\n" % (op.filename, op.ref))
  out.write("# Std: strand of reads\n")
  out.write("Std\tLength\tOccurences \n")
  for std in tab:
    for i in tab[std]:
      # write Length in 1-base offset
      out.write("%s\t%d\t%d\n" % (std, i+1, tab[std][i]))


def dmgFreqIsLow(folder):
  """ Returns true if the damage frequencies are too low to allow
  Bayesian estimation of DNA damages, i.e < 1% at first position.
  """
  total = 0.0
  for filename in ("5pCtoT_freq.txt", "5pGtoA_freq.txt"):
    with open(os.path.join(folder, filename)) as handle:
      for line in handle:
        freq = line.strip().split('\t')
        if freq[0] == "1":
          total += float(freq[1])
          break
      else:
        print("Error: Could not find pos = 1 in table '%s', bayesian computation cannot be performed" \
              % filename)
        return True

  if total < 0.01:
    print("Warning: DNA damage levels are too low, bayesian computation cannot be performed (%f < 0.01)\n" % total)
    return True

  return False
