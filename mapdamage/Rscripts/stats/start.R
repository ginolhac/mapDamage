logLikAllOptimize  <- function(x){
    #Optim wrapper for the log likelihood function
    #1 Theta
    #2 DeltaD
    #3 DeltaS
    #4 Lambda
    #5 LambdaRight
    #6 Lambda disp
    if (sum(x[1:5]>1)>0 || sum(x[1:6]<0)>0){
        return(Inf)
    }
    theta_mat  <- getTheta(x[1])
    leftLaVec  <- seqProbVecLambda(x[4],x[6],cu_pa$m,cu_pa$forward_only,cu_pa$reverse_only)#This is ugly using global variables ....
    rightLaVec  <- seqProbVecLambda(x[5],x[6],cu_pa$m,cu_pa$forward_only,cu_pa$reverse_only)
    if (cu_pa$forward_only){
        #Only using the forward end
        laVec <- leftLaVec
    }else if (cu_pa$reverse_only){
        #Only using the backward end
        laVec <- rightLaVec
    }else {
        #Both ends
        laVec <- c(leftLaVec[1:(cu_pa$m/2)],rightLaVec[(cu_pa$m/2+1):cu_pa$m])
    }
    return(-logLikAll(cu_pa$dat,theta_mat,x[2],x[3],laVec,cu_pa$nuVec,cu_pa$m))
}


logLikAllOptimizeWOLambdaDisp <- function(x){
    #1 Theta
    #2 DeltaD
    #3 DeltaS
    #4 Lambda
    #5 LambdaRight
    return(logLikAllOptimize(c(x[1],x[2],x[3],x[4],x[5],cu_pa$LambdaDisp)))
}

logLikAllOptimizeWOSamOverhang <- function(x){
    #1 Theta
    #2 DeltaD
    #3 DeltaS
    #4 Lambda
    #5 LambdaDisp
    return(logLikAllOptimize(c(x[1],x[2],x[3],x[4],x[4],x[5])))
}

logLikAllOptimizeWOLambdaDispAndSamOverhang <- function(x){
    #1 Theta
    #2 DeltaD
    #3 DeltaS
    #4 Lambda
    return(logLikAllOptimize(c(x[1],x[2],x[3],x[4],x[4],cu_pa$LambdaDisp)))
}


gridSearch <- function(cp,iter){
    if (cp$nuSamples!=0){
        print("Can't use the grid search with simulated nu vector")
        stop()
    }
    cat("Starting grid search, starting from random values\n")
    minVal <- Inf
    minParams <- list()
    control  <- list(maxit=5000)
    for (i in 1:iter){
        if (cp$fix_disp && cp$same_overhangs){
            #Fixing the dispersion and overhangs are the same on both sides ...
            out <- optim(c(runif(2,max=0.05),runif(2)),logLikAllOptimizeWOLambdaDispAndSamOverhang,method="Nelder-Mead",control=control)
        }else if (cp$fix_disp && !cp$same_overhangs) {
            #Fixing the dispersion and overhangs are different ...
            out <- optim(c(runif(2,max=0.05),runif(3)),logLikAllOptimizeWOLambdaDisp,method="Nelder-Mead",control=control)
        }  else if (!cp$fix_disp && cp$same_overhangs){
            #Variable dispersion and over hangs are same 
            out <- optim(c(runif(2,max=0.05),runif(2),sample(c(0.5,1,2,3,4,50,100,150,400),1)),logLikAllOptimizeWOSamOverhang,method="Nelder-Mead",control=control)
        }  else if (!cp$fix_disp && !cp$same_overhangs){
            #Variable dispersion and over hangs are different (Most general case)
            out <- optim(c(runif(2,max=0.05),runif(3),sample(c(0.5,1,2,3,4,50,100,150,400),1)),logLikAllOptimize,method="Nelder-Mead",control=control)
        }else {
            print("There is something wrong with the settings")
            stop()
        }
        if (out$value<minVal){
            minParams <- out
            minVal <- out$value
        }
        cat("Iteration\t",i,"\t",-minVal,"\n")
    }
    cp$Theta <- minParams$par[1]
    cp$ThetaMat <- getTheta(cp$Theta)
    cp$DeltaD <- minParams$par[2]
    cp$DeltaS <- minParams$par[3]
    cp$Lambda <- minParams$par[4]
    cp$laVec <- seqProbVecLambda(cp$Lambda,cp$LambdaDisp,nrow(cp$dat),cp$forward_only,cp$reverse_only)
    cp$old_lik <- - minVal
    #Since optim function returns the parameters in different orders it makes 
    #little bit complicated.
    if (cp$fix_disp && cp$same_overhangs){
        #Fixing the dispersion and overhangs are the same on both sides ...
    }else if (cp$fix_disp && !cp$same_overhangs) {
        #Fixing the dispersion and overhangs are different ...
        cp$LambdaRight <- minParams$par[5]
        cp$laVecRight <- seqProbVecLambda(cp$LambdaRight,cp$LambdaDisp,nrow(cp$dat),cp$forward_only,cp$reverse_only)
    } else if (!cp$fix_disp && cp$same_overhangs){
        #Variable dispersion and over hangs are same 
        cp$LambdaDisp <- minParams$par[5]
    } else if (!cp$fix_disp && !cp$same_overhangs){
        #Variable dispersion and over hangs are different (Most general case)
        cp$LambdaRight <- minParams$par[5]
        cp$LambdaDisp <- minParams$par[6]
        cp$laVecRight <- seqProbVecLambda(cp$LambdaRight,cp$LambdaDisp,nrow(cp$dat),cp$forward_only,cp$reverse_only)
    }else {
        print("There is something wrong with the settings")
        stop()
    }
    return(cp)
}
