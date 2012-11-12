#!/usr/bin/env python

import os
import subprocess
from mapdamage.version import __version__

def plot(op):

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

#os.path.dirname(__file__)
  script = os.path.join(os.path.dirname(__file__), "Rscripts/mapDamage.R")
  call = ["Rscript", script, fcomp, title, op.refplot, fmut, op.readplot, op.ymax, op.folder, op.title, __version__]
  code = subprocess.call(map(str, call))
  #print " ".join(map(str, call))

  if code == 0:
    print("pdf %s generated using R" % title)
    return 0
  else:
    print("Error: plotting with R failed")
    return 0
      


