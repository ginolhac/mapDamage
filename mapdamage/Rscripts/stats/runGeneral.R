#!/usr/bin/Rscript
source("../common.r")

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

# Number of random starting points for the grid search
grid_iter <- getArgument("GRID_ITER", as.integer)
# Burn in period
burn_in <- getArgument("BURN_IN", as.integer)
# Adjust proposal variance parameters iterations
adjust_iter <- getArgument("ADJUST_ITER", as.integer)
# Iterations
iterations <- getArgument("ITERATIONS", as.integer)

# Taking only the 5' end of the seqs 
forward_only <- getArgument("FORWARD_ONLY", as.logical)
# Taking only the 3' end of the seqs 
reverse_only <- getArgument("REVERSE_ONLY", as.logical)
# Geom instead of neg Bin 
fix_disp <- getArgument("FIX_DISP", as.logical)
# The overhangs are the same on both sides
same_overhangs <- getArgument("SAME_OVERHANGS", as.logical)

# Estimate the nu vector using the Briggs model (Should be similiar ammount to the number of sequences use this on your own risk....)
nu_samples <- getArgument("NU_SAMPLES", as.integer)
# Set 1 at 5' end and 0 at 3' end or else estimates it with GAM
fix_nu <- getArgument("FIX_NU", as.logical)
# Single stranded protocol C>T at both sides
ds_protocol <- getArgument("DS_PROTOCOL", as.logical)
# How long sequence to use from each side 
sub_length <- getArgument("SUB_LENGTH", as.integer)

# Absolute path to the dataset
path_to_dat <- getArgument("PATH_TO_DAT")
# Base file name of the output
out_file_base  <- getArgument("OUT_FILE_BASE")
# These options control the volume of the output
verbose <- getArgument("VERBOSE", as.logical)
quiet <- getArgument("QUIET", as.logical)

# Fix the transition and transversion ratio and acgt frequencies are equal
jukes_cantor <- getArgument("JUKES_CANTOR", as.logical)
use_raw_nick_freq <- getArgument("USE_RAW_NICK_FREQ", as.logical)
use_bw_theme <- getArgument("USE_BW_THEME", as.logical)


#######################################################
#
#                   Run the program
#
#######################################################

source("main.R")
