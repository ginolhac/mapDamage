# Prior and proposal distributions for the parameter,
# they are all uninformative.

priorTheta <- function(x) {
  return(dnorm(x = x, mean = 1, sd = 500, log = TRUE))
}


priorRho <- function(x) {
  return(dnorm(x = x, mean = 1, sd = 500, log = TRUE))
}


priorDeltaD <- function(x) {
  if (x < 0 || x > 1) {
    return(-Inf)
  }
  return(dbeta(x = x, shape1 = 1, shape2 = 1, log = TRUE))
}


priorDeltaS <- function(x) {
  if (x < 0 || x > 1) {
    return(-Inf)
  }
  return(dbeta(x = x, shape1 = 1, shape2 = 1, log = TRUE))
}


priorLambda <- function(x) {
  if (x < 0 || x > 1) {
    return(-Inf)
  }
  return(dbeta(x = x, shape1 = 1, shape2 = 1, log = TRUE))
}


priorLambdaRight <- function(x) {
  if (x < 0 || x > 1) {
    return(-Inf)
  }
  return(dbeta(x = x, shape1 = 1, shape2 = 1, log = TRUE))
}


priorLambdaDisp <- function(x) {
  if (x < 0) {
    return(-Inf)
  }
  # 2 times since we truncate it at zero
  return(log(2) + dnorm(x = x, mean = 0, sd = 100, log = TRUE))
}


proposeTheta <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$Theta))
}


proposeRho <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$Rho))
}


proposeDeltaD <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$DeltaD))
}


proposeDeltaS <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$DeltaS))
}


proposeLambda <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$Lambda))
}


proposeLambdaRight <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$LambdaRight))
}


proposeLambdaDisp <- function(x = NA) {
  return(rnorm(1, mean = x, sd = proposeParameters$LambdaDisp))
}
