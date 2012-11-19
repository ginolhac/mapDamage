#!/usr/bin/env python

import os
import sys
import collections

import mapdamage
from mapdamage.version import __version__


def initializeMut(ref, lg):  
  tab = {}
  for contig in ref:
    tab_contig = tab[contig] = {}
    for end in ('5p','3p'):
      tab_end = tab_contig[end] = {}
      for std in ('+','-'):
        tab_std = tab_end[std] = {}
        for mut in mapdamage.seq.header:
          tab_std[mut] = dict.fromkeys(xrange(lg), 0)
  
  return tab


def printMut(mut, op, out):
  _print_freq_table(mut, mapdamage.seq.header, op, out, offset = 1)


def initializeComp(ref, around, lg):
  keys = {"3p" : range(-lg, 0) + range(1, around + 1),
          "5p" : range(-around, 0) + range(1, lg + 1)}

  tab = {}
  for contig in ref:
    tab_contig = tab[contig] = {}
    for end in ('5p','3p'):
      tab_end = tab_contig[end] = {}
      for std in ('+','-'):
        tab_std = tab_end[std] = {}
        for letters in mapdamage.seq.letters:
          tab_std[letters] = dict.fromkeys(keys[end], 0)

  return tab


def printComp(comp, op, out):
  columns = mapdamage.seq.letters + ("Total",)
  _print_freq_table(comp, columns, op, out)


def initializeLg():
  tab = {}
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


def _print_freq_table(table, columns, op, out, offset = 0):
  out.write("# table produced by mapDamage version %s\n" % __version__)
  out.write("# using mapped file %s and %s as reference file\n" % (op.filename, op.ref))
  out.write("# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads\n")
  out.write("Chr\tEnd\tStd\tPos\t%s\n" % ("\t".join(columns)))

  for (reference, ends) in sorted(table.iteritems()):
    for (end, strands) in sorted(ends.iteritems()):
      for (strand, subtable) in sorted(strands.iteritems()):
        subtable["Total"] = {}
        for index in sorted(subtable[columns[0]]):
          subtable["Total"][index] = sum(subtable[letter][index] for letter in mapdamage.seq.letters)

          out.write("%s\t%s\t%s\t%d" % (reference, end, strand, index + offset))
          for base in columns:
            out.write("\t%d" % subtable[base][index])
          out.write("\n")
