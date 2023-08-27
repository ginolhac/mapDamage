# Various useful functions

isNearly <- function(target, current) {
    # Returns true if target is (nearly) equal to current
    return(isTRUE(all.equal(target, current)))
}

getPmat <- function(tmu,tv_ti_ratio,acgt) {
    # Returns the evolutionary substitution matrix
    if (any(acgt >= 1) || any(acgt <= 0)) {
        abort("The ACGT frequencies must be in the range 0 to 1")
    } else if (!isNearly(sum(acgt), 1)) {
        abort("The ACGT frequencies do not sum to 1")
    } else if (tv_ti_ratio <= 0) {
        abort("The transversion and transtition ratio cannot go under 0")
    }

    # Returns the substitution probability matrix.
    if (identical(tv_ti_ratio, 1) && identical(acgt, c(0.25, 0.25, 0.25, 0.25))) {
        return(jukesCantorPmat(tmu))
    } else {
        Q <- qmatHKY85(tmu, tv_ti_ratio, acgt)
        r  <- eigen(Q)
        B <- r$vectors
        E  <- diag(exp(r$values))

        #      Q        The eigen vector change of basis
        #  M  ->   M
        #
        #  ^       ^
        #  |B^-1   |B
        #     Eig
        #  N  ->   N

        # Little trick to avoid numerical difficulties
        out <- solve(a=t(B), b=E %*% t(B))

        rownames(out) <- c("A","C","G","T")
        colnames(out) <- c("A","C","G","T")
        return(out)
    }
}

jukesCantorPmat2 <- function(tmu){
    # Using the Juke-Cantor model
    names <- c("A", "C", "G", "T")
    return(matrix(1 / 4 - exp(-tmu) / 4, nrow=4, ncol=4, dimnames=list(names, names)) + diag(exp(-tmu), 4, 4))
}

qmatHKY85 <- function(tmu, tv_ti, acgt) {
    # HKY85 model
    #               |  sum_1          pi_c * tv_ti   pi_g           pi_t * tv_ti  |
    #               |  pi_a * tv_ti   sum_2          pi_g * tv_ti   pi_t          |
    #  Q = tmu  *   |  pi_a           pi_c * tv_ti   sum_3          pi_t * tv_ti  |
    #               |  pi_a * tv_ti   pi_c           pi_g * tv_ti   sum_4         |
    #
    Qmat <- rbind(
        c(0, tv_ti, 1, tv_ti) * acgt,
        c(tv_ti, 0, tv_ti, 1) * acgt,
        c(1, tv_ti, 0, tv_ti) * acgt,
        c(tv_ti, 1, tv_ti, 0) * acgt)

    return(tmu * (Qmat - diag(rowSums(Qmat))))
}

metroDesc <- function(lpr, lol){
    # the logic in the Metropolis-Hastings step
    stopifnot(!is.na(lpr))
    stopifnot(!is.na(lol))

    return(log(runif(1)) < lpr - lol)
}

seqProbVecLambda <- function(lambda, lambda_disp, m, termini="both") {
    # Returns the position specific probability of being in an overhang
    pvals <- dnbinom(1:m - 1, prob=lambda, size=lambda_disp)
    psum <- (1 - cumsum(pvals)) / 2

    if (termini == "both") {
        return(c(psum[1:(m / 2)], rev(psum[1:(m / 2)])))
    } else if (termini == "5p") {
        return(psum)
    } else if (termini == "3p") {
        return(rev(psum))
    } else {
        abort("Invalid value for termini: '%s'", termini)
    }
}


cat("Compiling C++ extension. This may take a moment ..\n")
sourceCpp(code='
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>

#include <gsl/gsl_sf_gamma.h>

// [[Rcpp::export]]
double logLikFunOneBaseFast(
    const Rcpp::NumericVector Gen,
    const Rcpp::NumericVector S,
    const Rcpp::NumericVector theta,
    const double deltad,
    const double deltas,
    const Rcpp::NumericVector laVec,
    const Rcpp::NumericVector nuVec,
    const int m,
    const int lin
) {
    Rcpp::NumericVector pDam(4);
    double result = 0;

    for (int i = 0; i < laVec.size(); i++) {
        double la = laVec[i];
        double nu = nuVec[i];
        double pct = nu * (la * deltas + deltad * (1 - la));
        double pga = (1 - nu) * (la * deltas + deltad * (1 - la));

        pDam[0] = theta(lin - 1, 0) * 1 + theta(lin - 1, 2) * pga;
        pDam[1] = theta(lin - 1, 1) * (1 - pct);
        pDam[2] = theta(lin - 1, 2) * (1 - pga);
        pDam[3] = theta(lin - 1, 1) * pct + theta(lin - 1, 3) * 1;

        double p1 = gsl_sf_lnfact(Gen(i)) -
            gsl_sf_lnfact(S(i, 0)) -
            gsl_sf_lnfact(S(i, 1)) -
            gsl_sf_lnfact(S(i, 2)) -
            gsl_sf_lnfact(S(i, 3));

        double p2 = S(i, 0) * log(pDam[0]) +
            S(i, 1) * log(pDam[1]) +
            S(i, 2) * log(pDam[2]) +
            S(i, 3) * log(pDam[3]);

        result += p1 + p2;
    }
    return result;
}')
cat("Compilation of C++ extension completed\n")


logLikAll <- function(dat, Theta, deltad, deltas, laVec, nuVec, m) {
    # Calculates the logLikelihood for all the bases
    if (deltad < 0 || deltad > 1 || deltas < 0 || deltas > 1) {
        return(-Inf)
    }

    Asub <- dat[, "A.C"] + dat[, "A.G"] + dat[, "A.T"]
    ALL <- logLikFunOneBaseFast(dat[,"A"], cbind(dat[,"A"]-Asub,dat[,"A.C"],dat[,"A.G"],dat[,"A.T"]),Theta,deltad,deltas,laVec,nuVec,m,1)

    Csub <- dat[, "C.A"] + dat[, "C.G"] + dat[, "C.T"]
    CLL <- logLikFunOneBaseFast(dat[,"C"], cbind(dat[,"C.A"],dat[,"C"]-Csub,dat[,"C.G"],dat[,"C.T"]),Theta,deltad,deltas,laVec,nuVec,m,2)

    Gsub <- dat[, "G.A"] + dat[, "G.C"] + dat[, "G.T"]
    GLL <- logLikFunOneBaseFast(dat[,"G"], cbind(dat[,"G.A"],dat[,"G.C"],dat[,"G"]-Gsub,dat[,"G.T"]),Theta,deltad,deltas,laVec,nuVec,m,3)

    Tsub <- dat[, "T.A"] + dat[, "T.C"] + dat[, "T.G"]
    TLL <- logLikFunOneBaseFast(dat[,"T"], cbind(dat[,"T.A"],dat[,"T.C"],dat[,"T.G"],dat[,"T"]-Tsub),Theta,deltad,deltas,laVec,nuVec,m,4)

    return(ALL + CLL + GLL + TLL)
}


getParams <- function(cp){
    #Utility function nice to update the MCMC iterations matrix
    return(c(cp$Theta,cp$Rho,cp$DeltaD,cp$DeltaS,cp$Lambda,cp$LambdaRight,cp$LambdaDisp))
}

plotTrace<- function(dat,main,k=111){
    #Running median of the MCMC iterations
    plot(1:length(dat),dat,xlab="Iteration",ylab="",main=main,type="l")
}

plotEverything <- function(mcmcOut,hi=0,pl,thin=100){
    #Plots the MCMC traceplot in the form of a running median and
    #histogram of the MCMC iterations
    if (sum(c(cu_pa$same_overhangs==FALSE,
                    cu_pa$fix_disp==FALSE,
                    cu_pa$fix_ti_tv==FALSE))>1){
        #Check if I need to add a extra row
        a_extra_row <- 1
    }else {
        a_extra_row <- 0
    }
    par(mfrow=c(3,2+a_extra_row))
    if(hi){
        hist(mcmcOut$out[,"Theta"],main=expression(theta),xlab="",freq=FALSE)
        if (!mcmcOut$cu_pa$fix_ti_tv){
            hist(mcmcOut$out[,"Rho"],main=expression(rho),xlab="",freq=FALSE)
        }
        hist(mcmcOut$out[,"DeltaD"],main=expression(delta[d]),xlab="",freq=FALSE)
        hist(mcmcOut$out[,"DeltaS"],main=expression(delta[s]),xlab="",freq=FALSE)
        hist(mcmcOut$out[,"Lambda"],main=expression(lambda),xlab="",freq=FALSE)
        if (!mcmcOut$cu_pa$same_overhangs){
            hist(mcmcOut$out[,"LambdaRight"],main=expression(lambda[r]),xlab="",freq=FALSE)
        }
        if (!mcmcOut$cu_pa$fix_disp){
            hist(mcmcOut$out[,"LambdaDisp"],main=expression(sigma[lambda]),xlab="",freq=FALSE)
        }
        hist(mcmcOut$out[,"LogLik"],main="LogLik",xlab="",freq=FALSE)
    }else {
        plotTrace(mcmcOut$out[,"Theta"],main=expression(theta))
        if (!mcmcOut$cu_pa$fix_ti_tv){
            plotTrace(mcmcOut$out[,"Rho"],main=expression(theta))
        }
        plotTrace(mcmcOut$out[,"DeltaD"],main=expression(delta[d]))
        plotTrace(mcmcOut$out[,"DeltaS"],main=expression(delta[s]))
        plotTrace(mcmcOut$out[,"Lambda"],main=expression(lambda))
        if (!mcmcOut$cu_pa$same_overhangs){
            plotTrace(mcmcOut$out[,"LambdaRight"],main=expression(lambda[r]))
        }
        if (!mcmcOut$cu_pa$fix_disp){
            plotTrace(mcmcOut$out[,"LambdaDisp"],main=expression(sigma[lambda]))
        }
        plotTrace(mcmcOut$out[,"LogLik"],main="LogLik")
    }
    par(mfrow=c(1,1))
}

accRat <- function(da){
    #A rough measure of the acceptance ratio
    return(length(unique(da))/length(da))
}

adjustPropVar <- function(mcmc,propVar){
    #Adjust the proposal variance to get something near .22
    for (i in colnames(mcmc$out)){
        if (i=="LogLik"){
            next
        } else if (i=="LambdaRight" & mcmc$cu_pa$same_overhangs){
            next
        } else if (i=="LambdaDisp" & mcmc$cu_pa$fix_disp){
            next
        } else if (i=="Rho" & mcmc$cu_pa$fix_ti_tv){
            next
        }
        rat <- accRat(mcmc$out[,i])
        if (rat<0.1){
            propVar[[i]] <- propVar[[i]]/2
        } else if (rat>0.3) {
            propVar[[i]] <- propVar[[i]]*2
        }
    }
    return(propVar)
}

runGibbs <- function(cu_pa,iter){
    #Sampling over the conditional posterior distribution
    #for the parameters.
    esti <- matrix(nrow=iter,ncol=8)
    colnames(esti) <- c("Theta","Rho","DeltaD","DeltaS","Lambda","LambdaRight","LambdaDisp","LogLik")
    for (i in 1:iter){
        cu_pa<-updateTheta(cu_pa)
        if (!cu_pa$fix_ti_tv){
            #Fix the transition and transversion ratio
            cu_pa<-updateRho(cu_pa)
        }
        cu_pa<-updateDeltaD(cu_pa)
        cu_pa<-updateDeltaS(cu_pa)
        cu_pa<-updateLambda(cu_pa)
        if (!cu_pa$same_overhangs){
            #Not the same overhangs update lambda right
            cu_pa <- updateLambdaRight(cu_pa)
        }
        if (!cu_pa$fix_disp){
            #Allowing dispersion in the overhangs
            cu_pa<-updateLambdaDisp(cu_pa)
        }
        esti[i,c(1:7)] <- getParams(cu_pa)
        esti[i,"LogLik"] <- logLikAll(cu_pa$dat,cu_pa$ThetaMat,cu_pa$DeltaD,cu_pa$DeltaS,cu_pa$laVec,cu_pa$nuVec,cu_pa$m)
        if (! (i %% 1000) && cu_pa$verbose){
            cat("MCMC-Iter\t",i,"\t",esti[i,"LogLik"],"\n")
        }
    }
    return(list(out=esti,cu_pa=cu_pa))
}


simPredCheck <- function(da,output){
    #Simulates one draw from the posterior predictive distribution
    #and the probability of a C>T substitution because of a cytosine
    #demnation.
    bases <- da[,c("A","C","G","T")]
    #Constructing the lambda vector
    if (output$cu_pa$same_overhangs){
        laVec <- seqProbVecLambda(sample(output$out[,"Lambda"],1),
                                  sample(output$out[,"LambdaDisp"],1),
                                  output$cu_pa$m,
                                  output$cu_pa$termini)
    }else {
        laVecLeft <- seqProbVecLambda(sample(output$out[,"Lambda"],1),
                                      sample(output$out[,"LambdaDisp"],1),
                                      output$cu_pa$m)
        laVecRight <- seqProbVecLambda(sample(output$out[,"LambdaRight"],1),
                                       sample(output$out[,"LambdaDisp"],1),
                                       output$cu_pa$m)
        laVec <- c(laVecLeft[1:(output$cu_pa$m/2)],laVecRight[(output$cu_pa$m/2+1):output$cu_pa$m])
    }
    #Constructing the nu vector
    nuVec <- output$cu_pa$nuVec

    #Sample the other parameters
    des <- sample(output$out[,"DeltaS"],1)
    ded <- sample(output$out[,"DeltaD"],1)
    the <- sample(output$out[,"Theta"],1)
    rho <- sample(output$out[,"Rho"],1)
    pmat <- getPmat(the,rho,output$cu_pa$acgt)
    ptransCT <- pmat["C","T"]
    ptransCC <- pmat["C","C"]
    ptransGA <- pmat["G","A"]
    ptransGG <- pmat["G","G"]
    #
    coln <- c("A.C","A.G","A.T","C.A","C.G","C.T","G.A","G.C","G.T","T.A","T.C","T.G")
    subs <- matrix(NA,nrow=nrow(output$cu_pa$dat),ncol=4+length(coln))
    colnames(subs) <- c("A","C","G","T",coln)
    #
    damProb <- rep(NA,nrow(output$cu_pa$dat))
    damProbGA <- damProb
    for (i in 1:nrow(output$cu_pa$dat)){
        #Construct the site specific probabilities
        pct <- nuVec[i]*(laVec[i]*des+ded*(1-laVec[i]))
        pga <- (1-nuVec[i])*(laVec[i]*des+ded*(1-laVec[i]))
        pDamMat <- matrix(c(
                         1,0,0,0,
                         0,1-pct,0,pct,
                         pga,0,1-pga,0,
                         0,0,0,1
                         ),nrow=4,byrow=TRUE)
        ThetapDam <- pDamMat %*% pmat
        #Calculate the probability C.T due to cytosine demanation
        damProb[i] <- ptransCC*pct/(ptransCC*pct+ptransCT)
        #Do not forget the reverse complement
        damProbGA[i] <- ptransGG*pga/(ptransGG*pga+ptransGA)
        #Then draw from a multinomial distribution
        subs[i,c("A.C","A.G","A.T")] <- t(rmultinom(1,output$cu_pa$dat[i,"A"],ThetapDam[1,]))[-1]/output$cu_pa$dat[i,"A"]
        subs[i,c("C.A","C.G","C.T")] <- t(rmultinom(1,output$cu_pa$dat[i,"C"],ThetapDam[2,]))[-2]/output$cu_pa$dat[i,"C"]
        subs[i,c("G.A","G.C","G.T")] <- t(rmultinom(1,output$cu_pa$dat[i,"G"],ThetapDam[3,]))[-3]/output$cu_pa$dat[i,"G"]
        subs[i,c("T.A","T.C","T.G")] <- t(rmultinom(1,output$cu_pa$dat[i,"T"],ThetapDam[4,]))[-4]/output$cu_pa$dat[i,"T"]
    }
    return(list(subs=subs,damProb=damProb,damProbGA=damProbGA))
}

calcSampleStats <- function(da,X){
    #Summary statistics of the posterior distributions
    return(data.frame(
                      x=1:nrow(da),
                      pos=da[,"Pos"],
                      mea=apply(X,1,mean),
                      med=apply(X,1,median),
                      loCI=apply(X,1,quantile,c(0.025), na.rm = TRUE),
                      hiCI=apply(X,1,quantile,c(0.975), na.rm = TRUE)
                      ))
}

postPredCheck <- function(da,output,samples=10000){
    #Plots the 95% posterior predictive intervals with the data as lines.
    #Returns the site specific probability of a C>T or G>A substitution
    #because of a cytosine demnation.
    CTs <- matrix(NA,nrow=nrow(da),ncol=samples)
    GAs <- matrix(NA,nrow=nrow(da),ncol=samples)
    REs <- matrix(NA,nrow=nrow(da),ncol=samples)
    C2TProbs <- matrix(NA,nrow=nrow(da),ncol=samples)
    G2AProbs <- matrix(NA,nrow=nrow(da),ncol=samples)
    #Two indices here, 1-based and the relative one
    da <- cbind(1:(nrow(da)),da)
    colnames(da)[1] <- "oneBased"
    #Get the breaks depending on parameters for pretty plots
    if (output$cu_pa$termini == "both") {
        bres <- c(seq(from=1,to=floor(nrow(da)/2),by=2),rev(seq(from=nrow(dat),to=floor(nrow(da)/2)+1,by=-2)))
    } else if (output$cu_pa$termini == "5p") {
        bres <- seq(from=1,to=nrow(da),by=2)
    } else if (output$cu_pa$termini == "3p") {
        bres <- seq(from=nrow(dat),to=1,by=-2)
    }else {
        abort("There is something fishy with the options")
    }
    labs <- dat[bres,"Pos"]
    #Sample from the posterior predicitive distibution
    for (i in 1:samples){
        sam <- simPredCheck(da,output)
        C2TProbs[,i] <- sam$damProb
        G2AProbs[,i] <- sam$damProbGA
        CTs[,i] <- sam$subs[,"C.T"]
        GAs[,i] <- sam$subs[,"G.A"]
        REs[,i] <- apply(sam$subs[,c("A.C","A.G","A.T","C.A","C.G","G.C","G.T","T.A","T.C","T.G")],1,mean)
    }
    CTsStats <- calcSampleStats(da,CTs)
    GAsStats <- calcSampleStats(da,GAs)
    REsStats <- calcSampleStats(da,REs)
    #Plotting the posterior predictive intervals
    p <- ggplot()+
         geom_point(aes(x,mea,colour="C->T"),data=CTsStats)+
         geom_point(aes(x,mea,colour="G->A"),data=GAsStats)+
         geom_point(aes(x,mea,colour="Others"),data=REsStats)+
         geom_errorbar(aes(x=x,ymin=loCI,ymax=hiCI,color="C->T"),data=CTsStats)+
         geom_errorbar(aes(x=x,ymin=loCI,ymax=hiCI,color="G->A"),data=GAsStats)+
         geom_errorbar(aes(x=x,ymin=loCI,ymax=hiCI,color="Others"),data=REsStats)+
         geom_line(aes(oneBased,C.T/C),color="red",data=data.frame(da))+
         geom_line(aes(oneBased,G.A/G),color="green",data=data.frame(da))+
         geom_line(aes(oneBased,((A.C+A.G+A.T)/A+(C.A+C.G)/C+(G.C+G.T)/G+(T.A+T.C+T.G)/T)/10),color="blue",data=data.frame(da))+
         labs(y = "Substitution rate",
              x = "Relative position",
              colour = "Subs. type",
              title = "Posterior prediction intervals")+
         scale_x_continuous(breaks=bres,labels=labs)
    if (output$cu_pa$use_bw_theme){
        p <- p+theme_bw()
    }
    plot(p)
    #The correcting probabilities
    coProbs <- cbind(da[,"Pos"],apply(C2TProbs,1,mean),apply(G2AProbs,1,mean))
    colnames(coProbs) <- c("Position","C.T","G.A")
    return(coProbs)
}


writeMCMC <- function(out,filename){
    #Writes the posterior samples to a file
    parameters <- c("Theta","DeltaD","DeltaS","Lambda")
    if (!out$cu_pa$fix_ti_tv){
        parameters <- c(parameters,"Rho")
    }
    if (!out$cu_pa$same_overhangs){
        parameters <- c(parameters,"LambdaRight")
    }
    if (!out$cu_pa$fix_disp){
        parameters <- c(parameters,"LambdaDisp")
    }
    parameters <- c(parameters,"LogLik")
    write.csv(out$out[,parameters],paste(filename,".csv",sep=""))
    #Now calculate summary statistic of the posterior distributions
    mea <- apply(out$out[,parameters],2,mean)
    std <- apply(out$out[,parameters],2,sd)
    qua <- apply(out$out[,parameters],2,quantile,seq(from=0,to=1,by=.025))
    acc <- apply(out$out[,parameters],2,accRat)
    summStat <- rbind(mea,std,acc,qua)
    rownames(summStat)[1] <- "Mean"
    rownames(summStat)[2] <- "Std."
    rownames(summStat)[3] <- "Acceptance ratio"
    write.csv(summStat,paste(filename,"_summ_stat.csv",sep=""))
}
