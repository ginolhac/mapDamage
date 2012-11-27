#!/usr/bin/env python

import os
import subprocess
from subprocess import CalledProcessError, check_call
from mapdamage.version import __version__
import mapdamage

def construct_R_path(name):
    """Construct a path to the R script given the name"""
    return(os.path.join(mapdamage.__path__[0], "Rscripts", name))

def plot(opt):
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

  fmut = opt.folder+"/"+"misincorporation.txt"
  fcomp = opt.folder+"/"+"dnacomp.txt"
  title = opt.folder+"/"+"Fragmisincorporation_"+opt.title+".pdf"

  script = construct_R_path("mapDamage.R") 
  call = ["Rscript", script, fcomp, title, opt.refplot, fmut, opt.readplot, \
      opt.ymax, opt.folder, opt.title, __version__]
  code = subprocess.call(map(str, call))

  if code == 0 and not opt.quiet:
    print("pdf %s generated" % title)
    return 0
  else:
    print("Error: plotting with R failed")
    return 1


        
def opt_plots(opt):       
  """ optional plots for length distribution
  and cumulative C>T mutations, per strand """
  
  fmut = opt.folder+"/"+"misincorporation.txt"
  flength = opt.folder+"/"+"lgdistribution.txt"
  output = opt.folder+"/"+"Length_"+opt.title+".pdf"

  script = construct_R_path("lengths.R") 
  call = ["Rscript", script, flength, output, fmut, opt.length, \
      opt.title, __version__]
  code = subprocess.call(map(str, call))
  if code == 0 and not opt.quiet:
    print("additional pdf %s generated" % output)
    return 0
  else:
    print("Error: plotting with R failed")
    return 1


def check_R_one_lib(name):
    """Checks if a necessary R library is here."""
    with open(os.devnull, "w") as devnull:
        rpa = construct_R_path("stats/checkLibraries.R")
        return subprocess.call(["Rscript", rpa, "--args", name],
                               stdout = devnull,
                               stderr = devnull)


def check_R_lib():
    """
    Checks if the necessary R libraries are here, signal 
    otherwise
    """
    libs = ["inline", "ggplot2", "gam", "Rcpp", "RcppGSL"]
    missing_lib = []
    for lib in libs:
        # Check the libraries
        if check_R_one_lib(lib):
            #found a missing library
            missing_lib.append(lib)
    if len(missing_lib) > 1:
        # Grammar Nazi has arrived
        last_ele = missing_lib.pop()
        print ("Missing the following R libraries '" + "', '".join(missing_lib) \
            + "' and '" + last_ele + "'")
        return 1
    elif len(missing_lib) == 1:
        print ("Missing the following R library "+missing_lib[0])
        return 1
    else :
        # No missing libraries
        return 0


def run_stats(opt):
    """
    Runs the Bayesian estimation program, using the options o
    """
    arg = ["Rscript",                            \
         construct_R_path("stats/runGeneral.R"), \
         "--args",                               \
         opt.rand,                               \
         opt.burn,                               \
         opt.adjust,                             \
         opt.iter,                               \
         opt.forward,                            \
         opt.reverse,                            \
         opt.fix_disp,                           \
         opt.same_hangs,                         \
         0,                                      \
         opt.fix_nicks,                          \
         opt.double_stranded,                    \
         opt.seq_length,                         \
         opt.folder+"/",                         \
         construct_R_path("stats/"),             \
         opt.folder+"/Stats_out",                \
         int(opt.verbose),                       \
         int(opt.quiet)
         ]  
    arg = [str(i) for i in arg]
    try:
        check_call(arg)
    except CalledProcessError:
        print("\nThe Bayesian statistics program failed to finish\n")
        raise SystemError
    return 0


