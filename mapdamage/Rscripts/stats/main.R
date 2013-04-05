#The main work flow of the package
#Load the libraries
suppressMessages(library(inline))  #Already checked the libraries 
suppressMessages(library(ggplot2)) #thus ignoring any messages from them
suppressMessages(library(Rcpp))
suppressMessages(library(gam)) 
suppressMessages(library(RcppGSL))

#Miscellaneous functions 
source(paste(path_to_mapDamage_stats,"function.R",sep=""))

#Prior and proposal distributions for the parameters
source(paste(path_to_mapDamage_stats,"priorPropose.R",sep=""))

#Update functions for the gibbs sampler
source(paste(path_to_mapDamage_stats,"postConditonal.R",sep=""))

#functions for the grid search
source(paste(path_to_mapDamage_stats,"start.R",sep=""))

#functions for loading the data 
source(paste(path_to_mapDamage_stats,"data.R",sep=""))


#######################################################
#
#              Loading the data
#
#######################################################

#This used for the Briggs MC estimation
if (nu_samples!=0){
    start_vals$lengths <- getSeqLen(path_to_dat)
    start_vals$le  <- floor(weighted.mean(x=start_vals$lengths$Length,w=start_vals$lengths$Occurences))
}

if (forward_only && reverse_only){
    write("Cannot specify using only the 5' end and the 3' end which makes no sense",stderr())
    stop()
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
              Nu=start_vals$nu,
              m=nrow(dat),
              meanLength=NA,
              lengths=NA,
              mLe=NA,
              laVec=NA,
              nuVec=NA,
              iter=iterations,
              nuSamples=nu_samples,
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
              quiet=quiet
              )

if (nu_samples!=0){
    #Need the lengths for the Briggs MC estimation
    cu_pa$meanLength <- start_vals$le
    cu_pa$lengths  <- start_vals$lengths
    cu_pa$mLe  <- max(start_vals$lengths$Length)
}

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
        write("Cannot use different overhangs with only the 5' end",stderr())
        stop()
    } else if (cu_pa$reverse_only){
        write("Cannot use different overhangs with only the 3' end",stderr())
        stop()
    }
    cu_pa$laVecRight <- seqProbVecLambda(cu_pa$LambdaRight,cu_pa$LambdaDisp,nrow(cu_pa$dat),cu_pa$forward_only,cu_pa$reverse_only)
}

#######################################################
#
#              Setting the nu vector 
#
#######################################################

if (!cu_pa$ds_protocol & cu_pa$nuSamples!=0){
    #Using the single stranded protocol with nu samples which makes 
    #no sense
    write("Silly to use the MC estimation for nu_i if using the single stranded protocol", stderr())
    stop()
}else if (cu_pa$nuSamples!=0 & cu_pa$fix_nu) {
    #Using the fixed nu parameter with nu samples which makes no sense 
    write("Silly to use the MC estimation for nu_i and use fix_nu ", stderr())
    stop()
}else if (cu_pa$nuSamples!=0){
    #Using the nu vector
    #IF ds protocol we set nu_vector as 1
    cu_pa$nuVec <- seqProbVecNuWithLengths(cu_pa$Lambda,
                                       cu_pa$LambdaDisp,
                               cu_pa$Nu,nrow(cu_pa$dat),
                               sampleHJ(cu_pa$lengths$Length,size=cu_pa$nuSamples,prob=cu_pa$lengths$Occurences),
                                              cu_pa$mLe,
                                     cu_pa$forward_only,
                                        cu_pa$nuSamples,
                                      cu_pa$ds_protocol)
}else if (!cu_pa$ds_protocol){
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
            cu_pa$nuVec  <- predict(gam(te~s(1:(cu_pa$m))))
        }else {
            cu_pa$nuVec  <- c(predict(gam(te[1:(cu_pa$m/2)]~s(1:(cu_pa$m/2)))),
                              predict(gam(te[(cu_pa$m/2+1):cu_pa$m]~s(1:(cu_pa$m/2)))))
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

if (out_file_base!=""){
    cat("Writing and plotting to files\n")
    #Print everything to file
    writeMCMC(mcmcOut,paste(out_file_base,"_MCMC_iter",sep=""))
    #Trace plot
    pdf(paste(out_file_base,"_MCMC_trace.pdf",sep=""))
    plotEverything(mcmcOut,hi=0)
    dev.off()
    #histogram of the conditional distributions
    pdf(paste(out_file_base,"_MCMC_hist.pdf",sep=""))
    plotEverything(mcmcOut,hi=1)
    dev.off()
    #Posterior predictive plot 
    pdf(paste(out_file_base,"_MCMC_post_pred.pdf",sep=""))
    siteProb <- postPredCheck(dat,mcmcOut)
    dev.off()
    #Write correcting probabilities 
    write.csv(siteProb,paste(out_file_base,"_MCMC_correct_prob.csv",sep=""))
} else {
    cat("Plotting\n")
    plotEverything(mcmcOut,hi=0)
    X11()
    plotEverything(mcmcOut,hi=1)
    X11()
    postPredCheck(dat,mcmcOut)
}
