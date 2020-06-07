#The main work flow of the package
#Load the libraries
suppressPackageStartupMessages(library(ggplot2)) #thus ignoring any messages from them
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(gam))
suppressPackageStartupMessages(library(RcppGSL))

#Miscellaneous functions
source("function.R")

#Prior and proposal distributions for the parameters
source("priorPropose.R")

#Update functions for the gibbs sampler
source("postConditonal.R")

#functions for the grid search
source("start.R")

#functions for loading the data
source("data.R")


#######################################################
#
#              Loading the data
#
#######################################################

if (forward_only && reverse_only){
    abort("Cannot specify using only the 5' end and the 3' end which makes no sense")
}

fow_dat <-readMapDamData(path_to_dat)
rev_dat <-readMapDamData(path_to_dat,forward=0)
if (forward_only){
    #Taking only the forward part
    dat <- fow_dat[1:sub_length,]
}else if (reverse_only){
    #Taking only the reverse part
    dat <- rev_dat[sub_length:1,]
}else {
    #Using both ends
    dat <- joinFowAndRev(fow_dat,rev_dat,sub_length)
}

#Getting everything ready for the mutation model
if (jukes_cantor){
    acgt <- c(0.25,0.25,0.25,0.25)
    fix_ti_tv <- TRUE
}else {
    acgt <- readBaseFreqs(path_to_dat)
    fix_ti_tv <- FALSE
}

#A list that keeps everything between the iterations
cu_pa <- list(
              dat=dat,
              ThetaMat=NA,
              Theta=-log((-start_vals$ptrans+.25)*4),
              Rho=start_vals$rho,
              DeltaD=start_vals$deltad,
              DeltaS=start_vals$deltas,
              Lambda=start_vals$lambda,
              LambdaRight=start_vals$lambda_right,
              LambdaDisp=start_vals$lambda_disp,
              m=nrow(dat),
              meanLength=NA,
              lengths=NA,
              mLe=NA,
              laVec=NA,
              nuVec=NA,
              iter=iterations,
              fix_nu=fix_nu,
              ds_protocol=ds_protocol,
              forward_only=forward_only,
              reverse_only=reverse_only,
              fix_disp= fix_disp,
              same_overhangs= same_overhangs,
              fix_ti_tv = fix_ti_tv,
              acgt = acgt,
              old_lik=-Inf,
              verbose=verbose,
              quiet=quiet,
              use_raw_nick_freq = use_raw_nick_freq,
              use_bw_theme = use_bw_theme
              )

cu_pa$ThetaMat <- getPmat(cu_pa$Theta,cu_pa$Rho,cu_pa$acgt)

#######################################################
#
#              Setting the lambda vector
#
#######################################################

cu_pa$laVec <- seqProbVecLambda(cu_pa$Lambda,cu_pa$LambdaDisp,nrow(cu_pa$dat),cu_pa$forward_only,cu_pa$reverse_only)

if (!cu_pa$same_overhangs){
    #The overhangs are not the same
    if (cu_pa$forward_only){
        abort("Cannot use different overhangs with only the 5' end")
    } else if (cu_pa$reverse_only){
        abort("Cannot use different overhangs with only the 3' end")
    }
    cu_pa$laVecRight <- seqProbVecLambda(cu_pa$LambdaRight,cu_pa$LambdaDisp,nrow(cu_pa$dat),cu_pa$forward_only,cu_pa$reverse_only)
}

#######################################################
#
#              Setting the nu vector
#
#######################################################

if (!cu_pa$ds_protocol){
    #The single stranded protocol
    cu_pa$nuVec <- rep(1,nrow(dat))
}else if (fix_nu){
    #Ones at the 5' end and zeros at the 3' end
    if (cu_pa$forward_only){
        cu_pa$nuVec <- rep(1,nrow(dat))
    }else if (cu_pa$reverse_only){
        cu_pa$nuVec <- rep(0,nrow(dat))
    }else {
        cu_pa$nuVec <- c(rep(1,nrow(dat)/2),rep(0,nrow(dat)/2))
    }
}else {
    #This is for a non linear nick frequency, assumes the G>T and G>A are
    #mostly due to DNA damage patterns do not use for low damage datasets
    te<-(dat[,"C.T"]/dat[,"C"])/(dat[,"G.A"]/dat[,"G"]+dat[,"C.T"]/dat[,"C"])
    if (sum(is.na(te) )!=0 ){
        write("Warning, To few substitutions to assess the nick frequency, using constant nick frequency instead", stderr())
        if (cu_pa$forward_only){
            cu_pa$nuVec <- rep(1,nrow(dat))
        }else if (cu_pa$reverse_only){
            cu_pa$nuVec <- rep(0,nrow(dat))
        }else {
            cu_pa$nuVec <- c(rep(1,nrow(dat)/2),rep(0,nrow(dat)/2))
        }
    } else {
        #The substitutes seem to be okay estimate the nick frequency using GAM
        if (cu_pa$forward_only || cu_pa$reverse_only){
            if (cu_pa$use_raw_nick_freq){
                #Use the frequency
                cu_pa$nuVec <- te
            }else{
                #Use a smoother
                cu_pa$nuVec  <- predict(gam(te~s(1:(cu_pa$m))))
            }
        }else {
            if (cu_pa$use_raw_nick_freq){
                cu_pa$nuVec  <- c(te[1:(cu_pa$m/2)],te[(cu_pa$m/2+1):cu_pa$m])
            }else{
                cu_pa$nuVec  <- c(predict(gam(te[1:(cu_pa$m/2)]~s(1:(cu_pa$m/2)))),
                                  predict(gam(te[(cu_pa$m/2+1):cu_pa$m]~s(1:(cu_pa$m/2)))))
            }
        }
        #This shouldn't happen but for sanity check
        cu_pa$nuVec[cu_pa$nuVec>1] <- 1
        cu_pa$nuVec[cu_pa$nuVec<0] <- 0
    }
    rm(te)
}
#######################################################
#
#          Finding an "optimal" starting place
#
#######################################################

if (grid_iter!=0){
    #Start at random places and optimize the likelihood function
    cu_pa <- gridSearch(cu_pa,grid_iter)
}

#Calculate the log likelihood in the beginning
if (cu_pa$same_overhangs){
    te_laVec <- cu_pa$laVec
}else {
    te_laVec <- c(cu_pa$laVec[1:(cu_pa$m/2)],cu_pa$laVecRight[(cu_pa$m/2+1):cu_pa$m])
}
cu_pa$old_lik <- logLikAll(cu_pa$dat,
                      cu_pa$ThetaMat,
                        cu_pa$DeltaD,
                        cu_pa$DeltaS,
                         te_laVec,
                         cu_pa$nuVec,
                             cu_pa$m
                  )


if (adjust_iter==0){
    #Single burning period
    if (!cu_pa$quiet){
        cat("Single burn in period\n")
    }
    mcmcOut <- runGibbs(cu_pa,burn_in)
    cu_pa <- mcmcOut$cu_pa
}else {
    for (i in 1:adjust_iter){
        if (!cu_pa$quiet){
            cat("Adjusting the proposal variance iteration ",i," \n")
        }
        mcmcOut <- runGibbs(cu_pa,burn_in)
        cu_pa <- mcmcOut$cu_pa
        proposeParameters <- adjustPropVar(mcmcOut,proposeParameters)
    }
}

if (!cu_pa$quiet){
    cat("Done burning, starting the iterations\n")
}
mcmcOut <- runGibbs(cu_pa,iterations)
cu_pa <- mcmcOut$cu_pa

if (!cu_pa$quiet){
    cat("Done with the iterations, finishing up\n")
}


cat("Writing and plotting to files\n")
# Print everything to file
writeMCMC(mcmcOut, paste(path_to_dat, "Stats_out_MCMC_iter", sep=""))

# Trace plot
pdf(paste(path_to_dat, "Stats_out_MCMC_trace.pdf", sep=""))
plotEverything(mcmcOut,hi=0)
dev.off()

# histogram of the conditional distributions
pdf(paste(path_to_dat, "Stats_out_MCMC_hist.pdf", sep=""))
plotEverything(mcmcOut,hi=1)
dev.off()

# Posterior predictive plot
pdf(paste(path_to_dat, "Stats_out_MCMC_post_pred.pdf", sep=""))
siteProb <- postPredCheck(dat,mcmcOut)
dev.off()

# Write correcting probabilities
write.csv(siteProb, paste(path_to_dat, "Stats_out_MCMC_correct_prob.csv", sep=""))
