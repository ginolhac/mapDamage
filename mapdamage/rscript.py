#!/usr/bin/env python

import os
import subprocess
from subprocess import CalledProcessError,check_call
from mapdamage.version import __version__
import mapdamage

def constructRPath (name):
    """Construct a path to the R script given the name"""
    return(os.path.join(mapdamage.__path__[0], "Rscripts",name))

def plot(op):
  """
  Calls the R script to make the plots, takes in the arguments 
  from parameter processing 
  """
  # comp<-args[1]
  # pdfout<-args[2]
  # around<-as.numeric(args[3])
  # misincorp<-args[4]
  # lg<-as.numeric(args[5])
  # ymax<-as.numeric(args[6])
  # folder<-args[7]
  # title<-args[8]
  # version<-args[9]

  fmut = op.folder+"/"+"misincorporation.txt"
  fcomp = op.folder+"/"+"dnacomp.txt"
  title = op.folder+"/"+"Fragmisincorporation_"+op.title+".pdf"

  script = constructRPath("mapDamage.R") 
  call = ["Rscript", script, fcomp, title, op.refplot, fmut, op.readplot, op.ymax, op.folder, op.title, __version__]
  code = subprocess.call(map(str, call))
  #print " ".join(map(str, call))

  if code == 0:
    print("pdf %s generated using R" % title)
    return 0
  else:
    print("Error: plotting with R failed")
    return 1

def checkRLib():
    """
    Checks if the necessary R libraries are here, terminate 
    the program otherwise
    """
    try:
      rpa=constructRPath("stats/checkLibraries.R")
      s=check_call(["Rscript",rpa],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      return 0
    except CalledProcessError:
      #Can't find the script or missing the libraries
      return 1

def runStats(o):
    """
    Runs the Bayesian estimation program, using the options o
    """
    arg=["Rscript",                      \
         constructRPath("stats/runGeneral.R"), \
         "--args",                       \
         o.rand,                         \
         o.burn,                         \
         o.adjust,                       \
         o.iter,                         \
         o.forward,                      \
         o.fix_disp,                     \
         o.same_hangs,                   \
         0,                              \
         o.fix_nicks,                    \
         o.double_stranded,              \
         o.seq_length,                   \
         o.folder+"/",                   \
         constructRPath("stats/"),             \
         o.folder+"/Stats_out"]  
    arg = [str(i) for i in arg]
    try:
        check_call(arg)
    except CalledProcessError:
        print("\nThe Bayesian statistics program failed to finish\n")
        raise SystemError
    return 0


