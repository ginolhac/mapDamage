#!/usr/bin/env python

import mapdamage
import sys 


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
  
  if read.is_reverse: 
    fivep = read.aend
    threep = read.pos
  else:
    fivep = read.pos
    threep = read.aend
  return(fivep, threep)


def getAround(coord, chrom, reflengths, lg, ref):
  
  i = j = 0
  before = after = ""
  if min(coord) - lg > 0:
    i = min(coord) - lg
  if max(coord) + lg > reflengths[chrom]:
    j = reflengths[chrom]
  else:
    j = max(coord) + lg 

  before = ref.fetch(chrom, i, min(coord))
  after = ref.fetch(chrom, max(coord), j)
  
  return(before, after)


def align(cigarlist, seq, ref):

  ins = parseCigar(cigarlist, 1)
  lref=list(ref)
  for nb,idx in ins:
    lref[idx:idx] = ["-"]*nb 
  ref = "".join(lref)
  delet = parseCigar(cigarlist, 2)
  lread = list(seq)
  for nb,idx in delet:
    lread[idx:idx] = ["-"]*nb 
  seq = "".join(lread)
  return(seq, ref)


def getMis(read, seq, refseq, ref, length, tab, end):
  if read.is_reverse:
    std = '-'
  else:
    std = '+'
  for (i, nt) in enumerate(zip(seq, refseq)):
    if i > length-1:
      continue
    # print("pos %d compare %s and %s" % (i, nt[0], nt[1]))
    mut = nt[1]+">"+nt[0] # mutation such as ref>read
    if nt[1] in mapdamage.seq.letters:
      # record base composition in the reference, only A, C, G, T
      tab[ref][end][std]['Tot'][i] += 1
      tab[ref][end][std][nt[1]][i] += 1
      #print("+1 for %s pos %d std %s (%d)" % (nt[1], i, std, misincorp[ref][std][nt[1]][i] ))
    try:
      # discard identities
      tab[ref][end][std][mut][i] += 1
      #print("+1 for %s pos %d std %s (%d)" % (mut, i, std, misincorp[ref][std][mut][i] ))
    except:
      continue


def parseCigar(cigarlist, op):
  tlength = 0
  coordinate = []
  # for deletion, count matches, softclip and deletion 
  if op == 2:
      oplist = (0, 2, 4, 7, 8)
  # for insertion and softclip count matches, softclip and insertions
  else:
      oplist = (0, 1, 4, 7, 8)
  for operation,length in cigarlist:
    if operation == op:
        coordinate.append([length, tlength])
    if operation in oplist: 
        tlength+=length
  return coordinate


# to be simplified without ref as arguments and subtables

def recordSoftClipping(sclip, read, ref, tab, lg):
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
          tab[ref]['3p']['-']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab[ref]['5p']['-']['S'][i] += 1
      else: 
        for i in range(0,min(sclip[idx][0], lg)):
          tab[ref]['5p']['+']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab[ref]['3p']['+']['S'][i] += 1
    # for soft clip at right side
    else:
      if read.is_reverse:
        for i in range(0,min(sclip[idx][0], lg)):
          tab[ref]['5p']['-']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab[ref]['3p']['-']['S'][i] += 1
      else:
        for i in range(0,min(sclip[idx][0], lg)):
          tab[ref]['3p']['+']['S'][i] += 1
        for i in range(min(len(read.seq)-sclip[idx][0], lg), min(len(read.seq), lg)):
          tab[ref]['5p']['+']['S'][i] += 1






