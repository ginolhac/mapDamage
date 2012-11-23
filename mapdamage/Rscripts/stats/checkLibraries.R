#!/usr/bin/Rscript
#Simple script to check if the necessary packages are installed.
args <- commandArgs(TRUE)
whichProgram = args[2]

if (whichProgram=="inline"){
    library(inline)
} else if (whichProgram=="ggplot2"){
    library(ggplot2)
} else if (whichProgram=="Rcpp"){
    library(Rcpp)
} else if  (whichProgram=="gam"){
    library(gam) 
} else if  (whichProgram=="RcppGSL"){
    library(RcppGSL)
} else {
    stop(paste("Error in checkLibraries.R, no need to check this library",whichProgram))
}
