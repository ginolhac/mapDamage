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

def checkROneLib(name):
    """Checks if a necessary R library is here."""
    with open(os.devnull, "w") as devnull:
        rpa = constructRPath("stats/checkLibraries.R")
        return subprocess.call(["Rscript", rpa, "--args", name],
                               stdout = devnull,
                               stderr = devnull)

def checkRLib():
    """
    Checks if the necessary R libraries are here, signal 
    otherwise
    """
    libs = ["inline","ggplot2","gam","Rcpp","RcppGSL"]
    missing_lib = []
    for lib in libs:
        #Check the libraries
        if checkROneLib(lib):
            #found a missing library
            missing_lib.append(lib)
    if (len(missing_lib)>1):
        #Grammar Nazi has arrived
        last_ele = missing_lib.pop()
        print ("Missing the following R libraries '" + "', '".join(missing_lib) + "' and '" + last_ele + "'")
        return 1
    elif (len(missing_lib)==1):
        print ("Missing the following R library "+missing_lib[0])
        return 1
    else :
        #No missing libraries
        return 0


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


