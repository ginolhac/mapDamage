#!/usr/bin/env python

import mapdamage
import sys 
import string
import itertools

# from Martin Kircher, description of CIGAR operations
#O BAM Description
#M 0 alignment match (can be a sequence match or mismatch)
#I 1 insertion to the reference
#D 2 deletion from the reference
#N 3 skipped region from the reference
#S 4 soft clipping (clipped sequences present in SEQ)
#H 5 hard clipping (clipped sequences NOT present in SEQ)
#P 6 padding (silent deletion from padded reference)
#= 7 sequence match
#X 8 sequence mismatch


def getCoordinates(read):  

  """ return external coordinates of aligned read bases """

  fivep = read.aend if read.is_reverse else read.pos
  threep = read.pos if read.is_reverse else read.aend
  
  return fivep, threep


def getAround(coord, chrom, reflengths, lg, ref):
  """ return reference sequences before and after the read
  check for extremities and return what is available """
  coord_min = min(coord)
  coord_max = max(coord)

  pos_before = max(0, coord_min - lg)
  pos_after  = min(reflengths[chrom], coord_max + lg)

  # Uppercased, to be sure that we don't compare A,C,G,T and a,c,g,t
  before = ref.fetch(chrom, pos_before, coord_min).upper()
  after = ref.fetch(chrom, coord_max, pos_after).upper()

  return before, after


def align(cigarlist, seq, ref):
  """ insert gaps according to the cigar string 
  deletion: gaps to be inserted into read sequences, 
  insertions: gaps to be inserted into reference sequence """  
  lref = list(ref)
  for nb, idx in parseCigar(cigarlist, 1):
    lref[idx:idx] = ["-"] * nb

  lread = list(seq)
  for nb, idx in parseCigar(cigarlist, 2):
    lread[idx:idx] = ["-"] * nb

  return "".join(lread), "".join(lref)


def getMis(read, seq, refseq, ref, length, tab, end):
  """ count mismatches using aligned reference and read,
  must be redone since N in reference were randomly replaced by any bases """
  std = '-' if read.is_reverse else '+'
  subtable = tab[ref][end][std]

  for (i, nt_seq, nt_ref) in itertools.izip(xrange(length), seq, refseq):
    if nt_ref in subtable:
      # record base composition in the reference, only A, C, G, T
      subtable[nt_ref][i] += 1

    # Most ref/seq pairs will be identical
    if (nt_ref != nt_seq):
      mut = "%s>%s" % (nt_ref, nt_seq)
      if mut in subtable:
        subtable[mut][i] += 1


def parseCigar(cigarlist, op):

  """ for a specific operation (mismach, match, insertion, deletion... see above
  return occurences and index in the alignment """
  tlength = 0
  coordinate = []
  # count matches, indels and mismatches
  oplist = (0, 1, 2, 7, 8)
  for operation,length in cigarlist:
    if operation == op:
        coordinate.append([length, tlength])
    if operation in oplist: 
        tlength+=length
  return coordinate


# to be simplified using subtables

def recordSoftClipping(sclip, read, tab, lg):
  #def update_table(left, right):
  #      subtable = tab[ref][left]['-']['S']
  #      for i in range(0,min(sclip[idx][0], lg)):
  #        subtable[i] += 1
  #      for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
  #        tab[ref][right]['-']['S'][i] += 1

  for idx in range(0,len(sclip)):
    # for soft at left side
    if idx == 0:
      if read.is_reverse:
        for i in range(0,min(sclip[idx][0], lg)):
          tab['3p']['-']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab['5p']['-']['S'][i] += 1
      else: 
        for i in range(0,min(sclip[idx][0], lg)):
          tab['5p']['+']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab['3p']['+']['S'][i] += 1
    # for soft clip at right side
    else:
      if read.is_reverse:
        for i in range(0,min(sclip[idx][0], lg)):
          tab['5p']['-']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab['3p']['-']['S'][i] += 1
      else:
        for i in range(0,min(sclip[idx][0], lg)):
          tab['3p']['+']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab['5p']['+']['S'][i] += 1






