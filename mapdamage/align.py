#!/usr/bin/env python

import itertools
import mapdamage

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


def get_coordinates(read):  
  """ return external coordinates of aligned read bases """
  fivep = read.aend if read.is_reverse else read.pos
  threep = read.pos if read.is_reverse else read.aend
  
  return fivep, threep


def get_around(coord, chrom, reflengths, length, ref):
  """ return reference sequences before and after the read
  check for extremities and return what is available """
  coord_min = min(coord)
  coord_max = max(coord)

  pos_before = max(0, coord_min - length)
  pos_after  = min(reflengths[chrom], coord_max + length)

  # Uppercased, to be sure that we don't compare A,C,G,T and a,c,g,t
  before = ref.fetch(chrom, pos_before, coord_min).upper()
  after = ref.fetch(chrom, coord_max, pos_after).upper()

  return before, after


def align(cigarlist, seq, ref):
  """ insert gaps according to the cigar string 
  deletion: gaps to be inserted into read sequences, 
  insertions: gaps to be inserted into reference sequence """  
  lref = list(ref)
  for nbr, idx in parse_cigar(cigarlist, 1):
    lref[idx:idx] = ["-"] * nbr

  lread = list(seq)
  for nbr, idx in parse_cigar(cigarlist, 2):
    lread[idx:idx] = ["-"] * nbr

  return "".join(lread), "".join(lref)


def align_with_qual(cigarlist, seq, qual, threshold, ref):
  """ insert gaps according to the cigar string 
  deletion: gaps to be inserted into read sequences and qualities, 
  insertions: gaps to be inserted into reference sequence """  
  lref = list(ref)
  for nbr, idx in parse_cigar(cigarlist, 1):
    lref[idx:idx] = ["-"] * nbr

  lread = list(seq)
  lqual = list(qual)
  for nbr, idx in parse_cigar(cigarlist, 2):
    lread[idx:idx] = ["-"] * nbr
    lqual[idx:idx] = ["-"] * nbr

  for idx, score in enumerate(lqual):
    # mask read and reference nucleotides (so not gaps) if below desired threshold
    if (ord(score)-33) < threshold and lread[idx] != "-":
      lread[idx] = "N"
      lref[idx]  = "N"

  return "".join(lread), "".join(lqual), "".join(lref)


def get_mis(read, seq, refseq, ref, length, tab, end):
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
      # discard masked nucleotides
      if mut in subtable:
        subtable[mut][i] += 1


def parse_cigar(cigarlist, ope):
  """ for a specific operation (mismach, match, insertion, deletion... see above)
  return occurences and index in the alignment """
  tlength = 0
  coordinate = []
  # count matches, indels and mismatches
  oplist = (0, 1, 2, 7, 8)
  for operation, length in cigarlist:
    if operation == ope:
        coordinate.append([length, tlength])
    if operation in oplist: 
        tlength += length
  return coordinate


def record_soft_clipping(read, tab, length):
  """ record soft clipped bases at extremities """
  def update_table(end, std, bases):
    for i in range(0, min(bases, length)):
      tab[end][std]['S'][i] += 1

  strand = '-' if read.is_reverse else '+'
  for (nbases, idx) in mapdamage.align.parse_cigar(read.cigar, 4):
    if idx == 0:
      # Soft-clipping at the left side of the alignment
      end = '3p' if read.is_reverse else '5p'
    else:
      # Soft-clipping at the right side of the alignment
      end = '5p' if read.is_reverse else '3p'

    update_table(end, strand, nbases)





