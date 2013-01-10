#! /usr/bin/Rscript
#Parses the command line arguments and calls the main program

rm(list=ls())
graphics.off()

argsList<-commandArgs(TRUE) #Assuming the user is not using this script directly
argsList <- argsList[-1]  #Skip --args


#######################################################
#
#                  Proposal parameters 
#
#######################################################

proposeParameters <- list( 
                          Theta=0.0003,
                          Rho=0.001,
                          DeltaD=0.001,
                          DeltaS=0.009,
                          Lambda=0.008,
                          LambdaRight=0.008,
                          LambdaDisp=0.015,
                          Nu=0.001
                          )

#######################################################
#
#                  Initial values
#
#######################################################


start_vals <- list(
                   ptrans = 0.00396/3,
                   rho = 1,
                   deltad = 0.0285,
                   deltas = 0.269,
                   lambda = 0.27,
                   lambda_right = 0.27,
                   lambda_disp = 1,
                   len = 60,
                   nu = 0.0645
                   )


#######################################################
#
#                   Various parameters
#
#######################################################

grid_iter <- as.integer(argsList[1])                 # Number of random starting points for the grid search
burn_in <- as.integer(argsList[2])                   # Burn in period
adjust_iter <- as.integer(argsList[3])               # Adjust proposal variance parameters iterations
iterations <- as.integer(argsList[4])                # Iterations

forward_only <- as.logical(as.numeric(argsList[5]))  # Taking only the 5' end of the seqs 
reverse_only <- as.logical(as.numeric(argsList[6]))  # Taking only the 3' end of the seqs 
fix_disp <- as.logical(as.numeric(argsList[7]))      # Geom instead of neg Bin 
same_overhangs <- as.logical(as.numeric(argsList[8]))# The overhangs are the same on both sides

nu_samples <- as.integer(argsList[9])                # Estimate the nu vector using the Briggs model (Should be similiar ammount to the number of sequences use this on your own risk....)
fix_nu <- as.logical(as.numeric(argsList[10]))       # Set 1 at 5' end and 0 at 3' end or else estimates it with GAM
ds_protocol <- as.logical(as.numeric(argsList[11]))  # Single stranded protocol C>T at both sides
sub_length <- as.integer(argsList[12])               # How long sequence to use from each side 

path_to_dat <- argsList[13]                          # Absolute path to the dataset
path_to_mapDamage_stats <- argsList[14]              # Absolute path to the mapDamage-stats folder
out_file_base  <- argsList[15]                       # Base file name of the output
verbose <- as.logical(as.numeric(argsList[16]))      # These options control the volume of the output
quiet <- as.logical(as.numeric(argsList[17]))

jukes_cantor <- as.logical(as.numeric(argsList[18])) # Fix the transition and transversion ratio and acgt frequencies are equal
path_to_acgt <- argsList[19] 


#######################################################
#
#                   Run the program
#
#######################################################

source(paste(path_to_mapDamage_stats,"main.R",sep=""))
