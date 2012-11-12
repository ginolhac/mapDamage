#!/usr/bin/env python

import string

# from Martin Kircher, to complement DNA
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv')

letters = ('A','C','G','T','Total')
mutations = ('G>A','C>T','A>G','T>C','A>C','A>T','C>G','C>A','T>G',\
        'T>A','G>C','G>T','A>-','T>-','C>-','G>-','->A','->T','->C','->G','S') 
header=letters+mutations

def writeFasta(read,ref, seq, refseq, start, end, before, after, fout):
  
  std = '-' if read.is_reverse else '+'

  # output coordinate in 1-based offset
  fout.write(">%s:%d-%d\n%s\n" % (ref, start-len(before), start, before))
  fout.write(">%s:%d-%d\n%s\n>%s_%s\n%s\n" % (ref, start+1, end+1, refseq, read.qname, std, seq))
  fout.write(">%s:%d-%d\n%s\n" % (ref, end+2, end+2+len(after), after))

def revcomp(seq):

  return seq.translate(table)[::-1]


