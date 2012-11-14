#!/usr/bin/env python

import string

# from Martin Kircher, to complement DNA
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv')

letters = ('A','C','G','T','Total')
mutations = ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G',\
        'T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') 
header=letters+mutations

def writeFasta(read, ref, seq, refseq, start, end, before, after, fout):
  std = '-' if read.is_reverse else '+'

  # output coordinate in 1-based offset
  if len(before) > 0:
    fout.write(">%s\n%s\n" % (ref, before))
  fout.write(">%s:%d-%d\n%s\n>%s_%s\n%s\n" % (ref, start+1, end+1, refseq, read.qname, std, seq))
  if len(after) > 0:
    fout.write(">%s\n%s\n" % (ref, after))


def revcomp(seq):
  """ return reverse complemented string """
  return seq.translate(table)[::-1]


def recordLg(read, coordinate, tab):
  """ record global length distribution
  don't record paired reads as they are normally not used for aDNA """

  std = '-' if read.is_reverse else '+'
  lg = (max(coordinate) - min(coordinate))
  if not read.is_paired:
    tab[std][lg] = tab[std][lg] + 1
    
  return tab


def fastaIndex(f):
  """ from a fasta index file, fai, return dictionary of references:lengths """
  
  fai = {}
  with open(f, 'r') as fh:
    for line in fh:
      ref = line.strip().split('\t')
      try:
        fai[ref[0]] = ref[1]
      except:
        sys.stderr.write("Error: fai file is not correct %s\n" % f)

  return fai


def checkRef(failengths, reflengths):
  """ check if all references in BAM/SAM are present in the fasta reference with same length """
  for ref in reflengths:
    lg = int(failengths.get(ref, 0))
    if lg != reflengths[ref]:
      print("\nreference %s is not correct in the fasta file" % ref)
      return None

  return True

