#Prior and proposal distributions for the parameter,
#they are all uninformative.

priorTheta <- function(x){
    return(dnorm(x=x,mean=1,sd=500,log=TRUE))
}

priorRho <- function(x){
    return(dnorm(x=x,mean=1,sd=500,log=TRUE))
}

priorDeltaD <- function(x){
    if (x<0 || x>1){
        return(-Inf)
    }
    return(dbeta(x=x,shape1=1,shape2=1,log=TRUE))
}

priorDeltaS <- function(x){
    if (x<0 || x>1){
        return(-Inf)
    }
    return(dbeta(x=x,shape1=1,shape2=1,log=TRUE))
}

priorLambda <- function(x){
    if (x<0 || x>1){
        return(-Inf)
    }
    return(dbeta(x=x,shape1=1,shape2=1,log=TRUE))
}

priorLambdaRight <- function(x){
    if (x<0 || x>1){
        return(-Inf)
    }
    return(dbeta(x=x,shape1=1,shape2=1,log=TRUE))
}

priorLambdaDisp <- function(x){
    if (x<0){
        return(-Inf)
    }
    #2 times since we truncate it at zero
    return(log(2)+dnorm(x=x,mean=0,sd=100,log=TRUE))
}

priorNu <- function(x){
    if (x<0 || x>1){
        return(-Inf)
    }
    return(dbeta(x=x,shape1=1,shape2=1,log=TRUE))
}

proposeTheta <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$Theta 
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeRho <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$Rho 
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeDeltaD <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$DeltaD
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)#Okay since norm is symmetric
    }
}

proposeDeltaS <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$DeltaS 
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeLambda <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$Lambda
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeLambdaRight <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$LambdaRight
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeLambdaDisp <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$LambdaDisp
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}

proposeNu <- function(x=NA,ra=NA){
    sh1 <- proposeParameters$Nu
    if (!is.na(ra)){
        return(rnorm(1,mean=x,sd=sh1))
    } else {
        return(0)
    }
}
