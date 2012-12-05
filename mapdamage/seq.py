#!/usr/bin/env python
import sys
import string

# from Martin Kircher, to complement DNA
TABLE = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb', \
    'ACGTKYWSRMBDHVacgtkywsrmbdhv')

LETTERS = ("A", "C", "G", "T")
MUTATIONS = ('G>A', 'C>T', 'A>G', 'T>C', 'A>C', 'A>T', 'C>G', 'C>A', 'T>G',
             'T>A', 'G>C', 'G>T', 'A>-', 'T>-', 'C>-', 'G>-', '->A', '->T', 
             '->C', '->G', 'S')
HEADER = LETTERS + ("Total", ) + MUTATIONS


def write_fasta(read, ref, seq, refseq, start, end, before, after, fout):
  std = '-' if read.is_reverse else '+'

  # output coordinate in 1-based offset
  if len(before) > 0:
    fout.write(">%s\n%s\n" % (ref, before))
  fout.write(">%s:%d-%d\n%s\n>%s_%s\n%s\n" % (ref, start+1, end+1, refseq, read.qname, std, seq))
  if len(after) > 0:
    fout.write(">%s\n%s\n" % (ref, after))


def revcomp(seq):
  """ return reverse complemented string """
  return seq.translate(TABLE)[::-1]


def record_lg(read, coordinate, tab):
  """ record global length distribution
  don't record paired reads as they are normally not used for aDNA """
  std = '-' if read.is_reverse else '+'
  
  length = (max(coordinate) - min(coordinate))
  if not read.is_paired:
    tab[std][length] = tab[std][length] + 1
    
  return tab


def read_fasta_index(filename):
  """ from a fasta index file, fai, return dictionary of references:lengths """
  def print_err(msg, filename, line):
    sys.stderr.write("Error: %s\n" % msg)
    sys.stderr.write("       Filename: %s\n" % filename)
    sys.stderr.write("       Line:     %s\n" % repr(line))
  
  fai = {}
  with open(filename, 'r') as handle:
    for line in handle:
      ref = line.split("\t")
      if len(ref) != 5:
        print_err("Line in fasta index contains wrong number of fields, found %i, expected 5:" \
            % len(ref), filename, line)
        return None

      try:
        fai[ref[0]] = int(ref[1])
      except ValueError:
        print_err("Column 2 in FASTA index did not contain a number, found '%s':" % ref[1], filename, line)
        return None

  return fai


def compare_sequence_dicts(fasta_dict, bam_dict):
  """Compares a FASTA and BAM sequence dictionary, and prints any differences.
  Returns true if all required sequences are found, false otherwise."""
  if fasta_dict == bam_dict:
    return True

  sys.stderr.write("Sequence dictionaries in FASTA/BAM files differ:\n")
  common = set(fasta_dict) & set(bam_dict)
  if not common:
    sys.stderr.write("FATAL ERROR: No sequences in common!\n")
    return False

  # Check that the lengths of required sequences match (fatal error)
  different = []
  for key in sorted(common):
    if fasta_dict[key] != bam_dict[key]:
      different.append((key, fasta_dict[key], bam_dict[key]))

  if different:
    sys.stderr.write("FATAL ERROR: Length of required sequences differ:\n")
    for values in different:
      sys.stderr.write("    - %s: %i bp vs %i bp\n" % values)

  # Check for sequences only found in the BAM file (fatal errors)
  bam_only = set(bam_dict) - common
  if bam_only:
    sys.stderr.write("FATAL ERROR: Sequences missing from FASTA dictionary:\n")
    for key in bam_only:
      sys.stderr.write("    - %s = %i bp\n" % (key, bam_dict[key]))

  # Check for sequences only found in the BAM file (fatal errors)
  fasta_only = set(fasta_dict) - common
  if fasta_only:
    sys.stderr.write("WARNING: FASTA file contains extra sequences:\n")
    for key in fasta_only:
      sys.stderr.write("    - %s = %i bp\n" % (key, fasta_dict[key]))

  sys.stderr.write("\n")

  return not (different or bam_only)
