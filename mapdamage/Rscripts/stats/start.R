#Runs likelihood optimization in the beginning to 
#start in a better place.

logLikAllOptimize  <- function(x,cp){
    #Optimizer wrapper for the log likelihood function
    Theta      <- x["Theta"]       
    DeltaD      <- x["DeltaD"]       
    DeltaS      <- x["DeltaS"]
    Lambda      <- x["Lambda"]
    LambdaRight <- x["LambdaRight"]
    LambdaDisp  <- x["LambdaDisp"]
    Rho         <- x["Rho"]
    if (sum(c(DeltaD,DeltaS,Lambda,LambdaRight)>1)>0 || sum(c(Theta,DeltaD,DeltaS,Lambda,LambdaRight,Rho)<0)>0){
        #This very convenient for lazy people, using a unbounded optimization method Nelder-Mead 
        #just returning Inf when we are outside the boundary.
        return(Inf)
    }
    if(cp$fix_ti_tv){
        theta_mat  <- getPmat(Theta,cp$Rho,cu_pa$acgt)
    }else {
        theta_mat  <- getPmat(Theta,Rho,cu_pa$acgt)
    }
    if (cp$fix_disp){
        disp_to_use =  cp$LambdaDisp
    }else {
        disp_to_use = LambdaDisp
    }
    leftLaVec  <- seqProbVecLambda(Lambda,disp_to_use,cu_pa$m,cu_pa$forward_only,cu_pa$reverse_only)
    rightLaVec  <- seqProbVecLambda(LambdaRight,disp_to_use,cu_pa$m,cu_pa$forward_only,cu_pa$reverse_only)
    if (cp$same_overhangs){
        rightLaVec <- leftLaVec
    }
    if (cp$forward_only){
        #Only using the forward end
        laVec <- leftLaVec
    }else if (cp$reverse_only){
        #Only using the backward end
        laVec <- rightLaVec
    }else {
        #Both ends
        laVec <- c(leftLaVec[1:(cu_pa$m/2)],rightLaVec[(cu_pa$m/2+1):cu_pa$m])
    }
    return(-logLikAll(cp$dat,theta_mat,DeltaD,DeltaS,laVec,cp$nuVec,cp$m))
}

gridSearch <- function(cp,iter){
    #Starts the Markov chain in local maxima for "quicker" burn-in period
    #No theoretical reasoning behind this but seems to work well in practice.
    if (cp$nuSamples!=0){
        write("Can't use the grid search with simulated nu vector",stderr())
        stop()
    }
    if(!cp$quiet){
        cat("Starting grid search, starting from random values\n")
    }
    minVal <- Inf
    minParams <- list()
    control  <- list(maxit=5000)
    for (i in 1:iter){
        out <- optim(
                     c(Theta=runif(1),
                       DeltaD=runif(1),
                       DeltaS=runif(1),
                       Lambda=runif(1),
                       LambdaRight=runif(1),
                       LambdaDisp=sample(c(0.5,1,2,3,4,50,100,150,400),1),
                       Rho=sample(c(0.5,.75,1,1.25,1.5),1)
                      ),
                     logLikAllOptimize,
                     NULL,
                     cp=cp,
                     method="Nelder-Mead",
                     control=control
                     )
        if (out$value<minVal){
            minParams <- out
            minVal <- out$value
        }
        if (cp$verbose){
            cat("Iteration\t",i,"\t",-minVal,"\n")
        }
    }
    #These parameters are always optimized
    cp$Theta <- minParams$par["Theta"]
    cp$DeltaD <- minParams$par["DeltaD"]
    cp$DeltaS <- minParams$par["DeltaS"]
    cp$Lambda <- minParams$par["Lambda"]
    
    #Only update the other ones if the user requested for them
    if (!cp$fix_ti_tv){
        cp$Rho  <- minParams$par["Rho"]
    }
    if (!cp$fix_disp){
        cp$LamdaDisp  <- minParams$par["LambdaDisp"]
    }
    if (!cp$same_overhangs){
        cp$LamdaRight  <- minParams$par["LambdaRight"]
    }

    #Now calculating matrix and vectors accompanying the optimal values
    cp$laVec <- seqProbVecLambda(cp$Lambda,cp$LambdaDisp,nrow(cp$dat),cp$forward_only,cp$reverse_only)
    cp$laVecRight <- seqProbVecLambda(cp$Lambda,cp$LambdaDisp,nrow(cp$dat),cp$forward_only,cp$reverse_only)
    cp$ThetaMat <- getPmat(cp$Theta,cp$Rho,cp$acgt)
    cp$old_lik <- - minVal

    return(cp)
}
