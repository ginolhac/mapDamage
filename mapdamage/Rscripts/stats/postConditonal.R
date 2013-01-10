#The posterior conditional function utilized by the 
# Gibbs sampler in function.R. They have all the 
#same form maybe they should be implemented in a 
#smarter way.

#The basic structure is the following 

#1. Get the old parameter 
#2. Propose a jump
#3. Accept it using the MH ratio
#4. Return the old or new value based on the MH ratio

updateTheta <- function(cp){
    old_lik  <- cp$old_lik+priorTheta(cp$Theta)
    theta_star <- proposeTheta(cp$Theta,1)
    if (theta_star<0 ){
        new_lik<- -Inf
    } else {
        theta_star_mat  <- getPmat(theta_star,cp$Rho,cp$acgt)
        new_lik_func  <-  logLikAll(cp$dat,theta_star_mat,cp$DeltaD,cp$DeltaS,cp$laVec,cp$nuVec,cp$m)
        new_lik <- new_lik_func+priorTheta(theta_star)
    }
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$Theta  <- theta_star
        cp$ThetaMat  <- theta_star_mat
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateRho <- function(cp){
    old_lik  <- cp$old_lik+priorRho(cp$Rho)
    rho_star <- proposeRho(cp$Rho,1)
    if (rho_star<=0 ){
        new_lik<- -Inf
    } else {
        rho_star_mat  <- getPmat(cp$Theta,rho_star,cp$acgt)
        new_lik_func  <-  logLikAll(cp$dat,rho_star_mat,cp$DeltaD,cp$DeltaS,cp$laVec,cp$nuVec,cp$m)
        new_lik <- new_lik_func+priorRho(rho_star)
    }
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$Rho  <- rho_star
        cp$ThetaMat  <- rho_star_mat
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateDeltaD <- function(cp){
    old_lik <- cp$old_lik+priorDeltaD(cp$DeltaD)
    deltad_star <- proposeDeltaD(cp$DeltaD,1)
    if (deltad_star<0 || deltad_star>1){
        new_lik <- -Inf
    }else {
        new_lik_func <- logLikAll(cp$dat,cp$ThetaMat,deltad_star,cp$DeltaS,cp$laVec,cp$nuVec,cp$m)
        new_lik <- new_lik_func+priorDeltaD(deltad_star)
    }
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$DeltaD  <- deltad_star
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateDeltaS <- function(cp){
    old_lik  <- cp$old_lik+priorDeltaS(cp$DeltaS)
    deltas_star <- proposeDeltaS(cp$DeltaS,1)
    if (deltas_star<0|| deltas_star>1){
        new_lik <- -Inf
    }else {
        new_lik_func  <- logLikAll(cp$dat,cp$ThetaMat,cp$DeltaD,deltas_star,cp$laVec,cp$nuVec,cp$m)
        new_lik  <- new_lik_func+priorDeltaS(deltas_star)
    }
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$DeltaS  <- deltas_star
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateLambda <- function(cp){
    old_lik  <- cp$old_lik+priorLambda(cp$Lambda)
    lambda_star <- proposeLambda(cp$Lambda,1)
    if (lambda_star<0 || lambda_star>1){
        return(cp)
    } 
    laVecStarLeft  <- seqProbVecLambda(lambda_star,cp$LambdaDisp,cp$m,cp$forward_only,cp$reverse_only)
    if (!cp$same_overhangs){
        #The left and right overhangs are not the same!
        laVecStar <- c(laVecStarLeft[1:(cp$m/2)],cp$laVecRight[(cp$m/2+1):cp$m])
    } else {
        #It left and right are the same
        laVecStar <- laVecStarLeft
    }
    new_lik_func  <- logLikAll(cp$dat,cp$ThetaMat,cp$DeltaD,cp$DeltaS,laVecStar,cp$nuVec,cp$m)
    new_lik  <- new_lik_func+priorLambda(lambda_star)
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$Lambda  <- lambda_star
        cp$laVec  <- laVecStar
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateLambdaRight <- function(cp){
    old_lik  <- cp$old_lik+priorLambdaRight(cp$LambdaRight)
    lambda_right_star <- proposeLambdaRight(cp$LambdaRight,1)
    if (lambda_right_star<0 || lambda_right_star>1){
        #Reject this right away
        return(cp)
    } 
    laVecStarRight  <- seqProbVecLambda(lambda_right_star,cp$LambdaDisp,cp$m,cp$forward_only,cp$reverse_only)
    if (!cp$same_overhangs){
        #The left and right overhangs are not the same!
        laVecStar <- c(cp$laVec[1:(cp$m/2)],laVecStarRight[(cp$m/2+1):cp$m])
    } else {
        #left and right are the same
        print("You shouldn't be calling this function if the overhangs are the same")
        stop()
    }
    if (cp$forward_only){
        print("You shouldn't be calling this function if you are only considering the forward part")
        stop()
    }
    if (cp$reverse_only){
        print("You shouldn't be calling this function if you are only considering the reverse part")
        stop()
    }
    new_lik_func  <- logLikAll(cp$dat,cp$ThetaMat,cp$DeltaD,cp$DeltaS,laVecStar,cp$nuVec,cp$m)
    new_lik  <- new_lik_func+priorLambdaRight(lambda_right_star)
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$LambdaRight  <- lambda_right_star
        cp$laVecRight  <- laVecStar
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateLambdaDisp <- function(cp){
    old_lik  <- cp$old_lik+priorLambdaDisp(cp$LambdaDisp)
    lambda_disp_star <- proposeLambdaDisp(cp$LambdaDisp,1)
    if (lambda_disp_star<0 ){
        return(cp)
    } 
    if (!cp$same_overhangs){
        leftLaVecStar  <- seqProbVecLambda(cp$Lambda,lambda_disp_star,cp$m,cp$forward_only,cp$reverse_only)
        rightLaVecStar  <- seqProbVecLambda(cp$LambdaRight,lambda_disp_star,cp$m,cp$forward_only,cp$reverse_only)
        laVecStar <- c(leftLaVecStar[1:(cp$m/2)],rightLaVecStar[(cp$m/2+1):cp$m])
    }else {
        laVecStar  <- seqProbVecLambda(cp$Lambda,lambda_disp_star,cp$m,cp$forward_only,cp$reverse_only)
    }
    new_lik_func  <- logLikAll(cp$dat,cp$ThetaMat,cp$DeltaD,cp$DeltaS,laVecStar,cp$nuVec,cp$m)
    new_lik  <- new_lik_func+priorLambdaDisp(lambda_disp_star)
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$LambdaDisp  <- lambda_disp_star
        cp$laVec  <- laVecStar
        cp$old_lik <- new_lik_func
    }
    return(cp)
}

updateNu <- function(cp){
    old_lik  <- cp$old_lik+priorNu(cp$Nu)
    nu_star <- proposeNu(cp$Nu,1)
    if (nu_star<0 || nu_star>1){
        return(cp)
    }
    nu_Vec_star <- seqProbVecNuWithLengths(cp$Lambda,cp$LambdaDisp,nu_star,cp$m,sampleHJ(cp$lengths$Length,size=cp$nuSamples,prob=cp$lengths$Occurences),cp$mLe,cp$forward_only,cu_pa$nuSamples,cu_pa$ds_protocol)
    new_lik_func  <- logLikAll(cp$dat,cp$ThetaMat,cp$DeltaD,cp$DeltaS,cp$laVec,nu_Vec_star,cp$m)
    new_lik  <- new_lik_func+priorNu(nu_star)
    if (metroDesc(new_lik,old_lik)) {
        #Accept
        cp$Nu  <- nu_star
        cp$nuVec  <- nu_Vec_star
        cp$old_lik <- new_lik_func
    }
    return(cp)
}
