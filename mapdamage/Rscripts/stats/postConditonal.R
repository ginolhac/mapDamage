# The posterior conditional function utilized by the
# Gibbs sampler in function.R. They have all the
# same form maybe they should be implemented in a
# smarter way.

# The basic structure is the following

# 1. Get the old parameter
# 2. Propose a jump
# 3. Accept it using the MH ratio
# 4. Return the old or new value based on the MH ratio

updateTheta <- function(cp) {
  theta_star <- proposeTheta(cp$Theta)
  if (theta_star < 0) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorTheta(cp$Theta)
  theta_star_mat <- getPmat(theta_star, cp$Rho, cp$acgt)
  new_lik_func <- logLikAll(cp$dat, theta_star_mat, cp$DeltaD, cp$DeltaS, cp$laVec, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorTheta(theta_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$Theta <- theta_star
    cp$ThetaMat <- theta_star_mat
    cp$old_lik <- new_lik_func
  }

  return(cp)
}


updateRho <- function(cp) {
  rho_star <- proposeRho(cp$Rho)
  if (rho_star <= 0) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorRho(cp$Rho)
  rho_star_mat <- getPmat(cp$Theta, rho_star, cp$acgt)
  new_lik_func <- logLikAll(cp$dat, rho_star_mat, cp$DeltaD, cp$DeltaS, cp$laVec, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorRho(rho_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$Rho <- rho_star
    cp$ThetaMat <- rho_star_mat
    cp$old_lik <- new_lik_func
  }

  return(cp)
}

updateDeltaD <- function(cp) {
  deltad_star <- proposeDeltaD(cp$DeltaD)
  if (deltad_star < 0 || deltad_star > 1) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorDeltaD(cp$DeltaD)
  new_lik_func <- logLikAll(cp$dat, cp$ThetaMat, deltad_star, cp$DeltaS, cp$laVec, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorDeltaD(deltad_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$DeltaD <- deltad_star
    cp$old_lik <- new_lik_func
  }

  return(cp)
}

updateDeltaS <- function(cp) {
  deltas_star <- proposeDeltaS(cp$DeltaS)
  if (deltas_star < 0 || deltas_star > 1) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorDeltaS(cp$DeltaS)
  new_lik_func <- logLikAll(cp$dat, cp$ThetaMat, cp$DeltaD, deltas_star, cp$laVec, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorDeltaS(deltas_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$DeltaS <- deltas_star
    cp$old_lik <- new_lik_func
  }

  return(cp)
}


updateLambda <- function(cp) {
  lambda_star <- proposeLambda(cp$Lambda)
  if (lambda_star < 0 || lambda_star > 1) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorLambda(cp$Lambda)
  laVecStarLeft <- seqProbVecLambda(lambda_star, cp$LambdaDisp, cp$m, cp$termini)
  if (cp$same_overhangs) {
    # It left and right are the same
    laVecStar <- laVecStarLeft
  } else {
    # The left and right overhangs are not the same
    laVecStar <- c(laVecStarLeft[1:(cp$m / 2)], cp$laVecRight[(cp$m / 2 + 1):cp$m])
  }

  new_lik_func <- logLikAll(cp$dat, cp$ThetaMat, cp$DeltaD, cp$DeltaS, laVecStar, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorLambda(lambda_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$Lambda <- lambda_star
    cp$laVec <- laVecStar
    cp$old_lik <- new_lik_func
  }

  return(cp)
}


updateLambdaRight <- function(cp) {
  stopifnot(!cp$same_overhangs, cp$termini == "both")

  lambda_right_star <- proposeLambdaRight(cp$LambdaRight)
  if (lambda_right_star < 0 || lambda_right_star > 1) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorLambdaRight(cp$LambdaRight)
  laVecStarRight <- seqProbVecLambda(lambda_right_star, cp$LambdaDisp, cp$m, cp$termini)
  # The left and right overhangs are not the same
  laVecStar <- c(cp$laVec[1:(cp$m / 2)], laVecStarRight[(cp$m / 2 + 1):cp$m])
  new_lik_func <- logLikAll(cp$dat, cp$ThetaMat, cp$DeltaD, cp$DeltaS, laVecStar, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorLambdaRight(lambda_right_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$LambdaRight <- lambda_right_star
    cp$laVecRight <- laVecStar
    cp$old_lik <- new_lik_func
  }

  return(cp)
}


updateLambdaDisp <- function(cp) {
  lambda_disp_star <- proposeLambdaDisp(cp$LambdaDisp)
  if (lambda_disp_star < 0) {
    return(cp)
  }

  old_lik <- cp$old_lik + priorLambdaDisp(cp$LambdaDisp)
  if (cp$same_overhangs) {
    laVecStar <- seqProbVecLambda(cp$Lambda, lambda_disp_star, cp$m, cp$termini)
  } else {
    leftLaVecStar <- seqProbVecLambda(cp$Lambda, lambda_disp_star, cp$m, cp$termini)
    rightLaVecStar <- seqProbVecLambda(cp$LambdaRight, lambda_disp_star, cp$m, cp$termini)
    laVecStar <- c(leftLaVecStar[1:(cp$m / 2)], rightLaVecStar[(cp$m / 2 + 1):cp$m])
  }
  new_lik_func <- logLikAll(cp$dat, cp$ThetaMat, cp$DeltaD, cp$DeltaS, laVecStar, cp$nuVec, cp$m)
  new_lik <- new_lik_func + priorLambdaDisp(lambda_disp_star)

  if (metroDesc(new_lik, old_lik)) {
    cp$LambdaDisp <- lambda_disp_star
    cp$laVec <- laVecStar
    cp$old_lik <- new_lik_func
  }

  return(cp)
}
