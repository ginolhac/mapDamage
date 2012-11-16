#!/usr/bin/env python

import mapdamage
import sys 
import string

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
  ins = parseCigar(cigarlist, 1)
  lref=list(ref) # to be sure that we don't compare A,C,G,T and a,c,g,t
  for nb,idx in ins:
    lref[idx:idx] = ["-"]*nb 
  ref = "".join(lref)
  delet = parseCigar(cigarlist, 2)
  lread = list(seq)
  for nb,idx in delet:
    lread[idx:idx] = ["-"]*nb 
  seq = "".join(lread)

  return seq, ref


def getMis(read, seq, refseq, ref, length, tab, end):

  """ count mismatches using aligned reference and read,
  must be redone since N in reference were randomly replaced by any bases """
  std = '-' if read.is_reverse else '+'
  
  for (i, nt) in enumerate(zip(seq, refseq)):
    if i > length-1:
      continue
    #print("pos %d compare %s and %s" % (i, nt[0], nt[1]))
    mut = nt[1]+">"+nt[0] # mutation such as ref>read
    if nt[1] in mapdamage.seq.letters:
      # record base composition in the reference, only A, C, G, T
      tab[ref][end][std]['Total'][i] += 1
      tab[ref][end][std][nt[1]][i] += 1
    try:
      # discard identities
      tab[ref][end][std][mut][i] += 1
    except:
      continue


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






