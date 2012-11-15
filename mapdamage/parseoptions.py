#!/usr/bin/env python
# -*- coding: ASCII -*-

from optparse import OptionParser,OptionGroup, SUPPRESS_HELP
import os
import sys

import mapdamage
from mapdamage.version import __version__
from mapdamage.rscript import checkRLib

def fileExist(filename):
  if os.path.exists(filename) and not os.path.isdir(filename):
    return True
  elif filename == "-":
    return True
  else:
    sys.stderr.write("Error: '%s' is not a valid file\n" % (filename))
    return None


def checkModule(mod):
  module = mod.keys()
  url = mod[module[0]]
  try:
    __import__(module[0])

  except ImportError, e:
    sys.stderr.write("Error: Could not import required module '%s':\n\t- %s\n" % (module[0],e))
    sys.stderr.write("       If module is not installed, please download from '%s'.\n" % (url,))
    sys.stderr.write("       A local install may be performed using the following command:\n")
    sys.stderr.write("       $ python setup.py install --user\n\n")
    return None

  return True


def whereis(program):
  for path in os.environ.get('PATH', '').split(':'):
    if os.path.exists(os.path.join(path, program)) and \
        not os.path.isdir(os.path.join(path, program)):
      return os.path.join(path, program)
  return None


def checkPyVersion():
  req_version = (2,6)
  cur_version = sys.version_info

  if cur_version >= req_version:
   return True
  else:
   sys.stderr.write("Your Python interpreter is too old."\
       "Please consider upgrading to at least %d.%d\n" % (req_version[0], req_version[1]))
   return None
    

def options(args):
  parser = OptionParser("%prog [options] -i BAMfile -r reference.fasta\n\nUse option -h or --help for help", version=__version__, \
          epilog="report bugs to aginolhac@snm.ku.dk, MSchubert@snm.ku.dk or jonsson.hakon@gmail.com")

  args = OptionGroup(parser, "Input files")
  args.add_option("-i", "--input", help="SAM/BAM file, must contain a valid header, use '-' for reading a BAM from stdin", \
        action="store", type="string", dest="filename")
  args.add_option("-r", "--reference", help="Reference file in FASTA format",\
        action="store", dest="ref")
 
  parser.add_option_group(args)
  group = OptionGroup(parser, "General options")
  group.add_option("-n", "--downsample", help = "Downsample to a randomly selected fraction of the reads (if 0 < DOWNSAMPLE < 1), or " \
                   "a fixed number of randomly selected reads (if DOWNSAMPLE >= 1). By default, no downsampling is performed.",
                   type = float, default = None)
  group.add_option("-l","--length",dest="length",help="read length, in nucleotides to consider [%default]",\
          type = int, default=70,action="store")
  group.add_option("-a","--around",dest="around",help="nucleotides to retrieve before/after reads [%default]",\
          type = int, default=10,action="store")
  group.add_option("-Q","--minbasequal",dest="minqual",help="minimun base quality Phred score considered [%default]",\
          type = int, default=0,action="store")
  #group.add_option("-p","--pipe",dest="pipe",help="Read BAM from a pipe, ie. the standard input",\
  #      default=False,action="store_true")
  group.add_option("-d", "--folder", help="folder name to store results [results_FILENAME]", \
        action="store", type="string", dest="folder")
  group.add_option("-f","--fasta",dest="fasta",help="Write alignments in a FASTA file",\
        default=False,action="store_true")
  group.add_option("--plotonly",dest="plotonly",help="Run only plotting from a valid result folder",\
        default=False,action="store_true")
  group.add_option("-q","--quiet",dest="quiet",help="Disable any output to stdout",\
        default=False,action="store_true")
  group.add_option("-v","--verbose",dest="verbose",help="Display progression information during parsing",\
        default=False,action="store_true")
  group.add_option("--noplot",dest="nor",help=SUPPRESS_HELP, default=False, action="store_true")
  parser.add_option_group(group)

  # options for plotting damage patterns
  group2 = OptionGroup(parser, "Options for graphics")
  group2.add_option("-y","--ymax",dest="ymax",\
          help="graphical y-axis limit for nucleotide misincorporation frequencies [%default]", type = float, default=0.3,action="store")
  group2.add_option("-m","--readplot",dest="readplot",\
          help="read length, in nucleotides, considered for plotting nucleotide misincorporations [%default]",\
          type = int, default=25,action="store")
  group2.add_option("-b","--refplot",dest="refplot",\
          help="the number of reference nucleotides to consider for ploting base composition in the region located upstream and downstream of every read [%default]",\
          type= int, default=10,action="store")
  group2.add_option("-t","--title",dest="title",\
          help="title used for both graph and filename [%default]",\
          type="string", default="plot",action="store")
  parser.add_option_group(group2)

  # Then the plethora of optional options for the statistical estimation ..
  group3 = OptionGroup(parser,"Options for the statistical estimation")
  group3.add_option("","--rand",dest="rand",\
          help="Number of random starting points for the likelihood optimization  [%default]", type = int, default=30,action="store")
  group3.add_option("","--burn",dest="burn",\
          help="Number of burnin iterations  [%default]", type = int, default=10000,action="store")
  group3.add_option("","--adjust",dest="adjust",\
          help="Number of adjust proposal variance parameters iterations  [%default]", type = int, default=10,action="store")
  group3.add_option("","--iter",dest="iter",\
          help="Number of final MCMC iterations  [%default]", type = int, default=50000,action="store")
  group3.add_option("","--forward",dest="forward",\
          help="Using only the 5' end of the seqs  [%default]", type = int, default=0,action="store")
  group3.add_option("","--fix_disp",dest="fix_disp",\
          help="Fix dispersion in the overhangs  [%default]", type = int, default=1,action="store")
  group3.add_option("","--same_hangs",dest="same_hangs",\
          help="The overhangs are the same on both sides  [%default]", type = int, default=1,action="store")
  group3.add_option("","--fix_nicks",dest="fix_nicks",\
          help="Fix the nick frequency vector nu else estimate it with GAM  [%default]", type = int, default=0,action="store")
  group3.add_option("","--double_stranded",dest="double_stranded",\
          help="Double stranded protocol [%default]", type = int, default=1,action="store")
  group3.add_option("","--seq_length",dest="seq_length",\
          help="How long sequence to use from each side [%default]", type = int, default=12,action="store")
  group3.add_option("--stats_only",dest="stats_only",help="Run only statistical estimation from a valid result folder",\
        default=False,action="store_true")
  group3.add_option("--nostats",dest="nos",help=SUPPRESS_HELP, default=False, action="store_true")

  parser.add_option_group(group3)

  #Parse the arguments
  (options, args) = parser.parse_args()
  
  # check python version
  if not checkPyVersion():
    return None

  # check general arguments
  if not (options.plotonly or options.stats_only) and not options.filename: 
    parser.error('SAM/BAM file not given')
  if not (options.plotonly or options.stats_only) and not options.ref:
    parser.error('Reference file not given')
  if not options.plotonly and not options.stats_only:
    if not fileExist(options.filename) or not fileExist(options.ref):
      return None
  if options.downsample is not None:
    if options.downsample <= 0:
      parser.error("-n/--downsample must be a positive value")
    elif options.downsample >= 1:
      options.downsample = int(options.downsample)
  
  if options.plotonly and not options.folder:
    parser.error('Folder not provided, required with --plotonly')
  if options.stats_only and not options.folder:
    parser.error('Folder not provided, required with --stats_only')

  # check options 
  if options.length < 0:
    parser.error('length (-l) must be a positive integrer')
  if options.around < 0:
    parser.error('around (-a) must be a positive integrer')
  if options.ymax <= 0 or options.ymax > 1:
    parser.error('ymax (-b) must be an integrer beetween 0 and 1')
  if options.readplot < 0:
    parser.error('readplot (-m) must be a positive integrer')
  if options.refplot < 0:
    parser.error('refplot (-b) must be a positive integrer')
  if options.refplot > options.around and not options.plotonly:
    parser.error('refplot (-b) must be inferior to around (-a)')
  if options.readplot > options.length:
    parser.error('readplot (-m) must be inferior to length (-l)')

  # check folder
  if not options.folder and options.filename:
    options.folder = "results_"+os.path.splitext(os.path.basename(options.filename))[0]
  if os.path.isdir(options.folder):
    if not options.quiet and not options.plotonly:
      print("Warning, %s already exists" % options.folder)
    if options.plotonly:
      if not fileExist(options.folder+"/dnacomp.txt") or not fileExist(options.folder+"/misincorporation.txt"):
        parser.error('folder %s is not a valid result folder' % options.folder)
  else:
    os.makedirs(options.folder, mode = 0700)
    if options.plotonly or options.stats_only:
      sys.stderr.write("Error, %s does not exist while plot/stats only was used\n" % options.folder)
      return None


  # check if the Rscript executable is present on the system
  if not whereis('Rscript'):
    print("Warning, Rscript is not in your PATH, plotting is disabled")
    options.nor = True
  
  if checkRLib():
    #Check for R libraries
    print("The Bayesian estimation is then disabled\n")
    options.nos = True

  return options

