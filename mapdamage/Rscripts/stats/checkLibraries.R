#!/usr/bin/Rscript
# Simple script to check if the necessary packages are installed.
for (name in commandArgs(TRUE)) {
    library(name, character.only=TRUE)
}
