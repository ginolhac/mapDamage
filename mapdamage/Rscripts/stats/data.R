sumTheChr <- function(da,co){
    #Sum up the columns for the different chromosomes
    return(tapply(da[[co]],da$Pos,sum))
}

readMapDamData <- function(folder,forward=1){
    #Reads in the data from mapdamage
    #Sums up the position counts for the different chromosomes
    fil <- "misincorporation.txt"
    raw_dat <- read.table(paste(folder,fil,sep=""),header=TRUE)
    if (forward==1){
        raw_dat <- raw_dat[raw_dat$End=="5p" & raw_dat$Std=="+",]
    } else {
        raw_dat <- raw_dat[raw_dat$End=="3p" & raw_dat$Std=="-",]
    }
    dat <- matrix(nrow=max(raw_dat$Pos),ncol=length(colnames(raw_dat))-3)
    colnames(dat) <- colnames(raw_dat)[c(-1,-2,-3)]
    dat[,"Pos"] <- seq(from=1,to=nrow(dat),by=1)
    for (i in colnames(raw_dat)[c(-1,-2,-3,-4)]){
        dat[,i] <- sumTheChr(raw_dat,i)
    }
    return(dat)
}

joinFowAndRev <- function(fo,re,nrPos){
    #Joins the 5' and 3' end with nrPos bases for each end
    te <- fo[1:nrPos,]
    te2 <- re[nrPos:1,]
    out <- rbind(te,te2)
    out[,"Pos"] <- seq(from=1,to=nrow(out),by=1)
    return(out)
}

getSeqLen <- function(pa){
    #Path to the mapDamage folder to get the length distribution
    le_dat_min <- read.table(paste(pa,"lengthDistribStrd-.txt",sep=""),header=TRUE)
    le_dat_plu <- read.table(paste(pa,"lengthDistribStrd+.txt",sep=""),header=TRUE)
    les <- list(Length=le_dat_min$Length,Occurences=le_dat_min$Occurences+le_dat_plu$Occurences)
    les$Length <- les$Length[les$Occurences!=0]
    les$Occurences <- les$Occurences[les$Occurences!=0]/sum(les$Occurences)
    return(les)
}



readMitoTestSet <- function(forward=TRUE){
    #This dataset was created in following fashion
    
    #head -n 2000 ../120124_700K_merged_bothaligned_Eqs_mito_noseed_noCloseHuman_kirdup.sam > first_2000.sam
    #../Damage192/damage_lik   -i first_2000.sam  -ref ../Equus_cab_mito3.fasta 2> /dev/stdout
    
    #got The following parameters using 
   #params: Delta   DeltaS  Lambda  JC      Nu
#mle:    0.025   0.851   0.417   0.0128  0.0719  -14906.46488    (converged) 
#14908.5	0.0251084 0.862156 0.422812 0.0128634 
#lowerbound:	0.0227	0.759	0.378	0.0122	0.0451
#upperbound:	0.0275	0.94	0.457	0.0135	0.132

    #Then the empirical distribution using 
    # ../Damage192/damage_lik -counts -empirical -includeSelf -i first_2000.sam  -ref ../Equus_cab_mito3.fasta
    
    dat <- read.table("horse_dataset/2000_reads_from_mito.txt",header=TRUE)
    
    if (forward){
        dat <- dat[1:24,]
    } else {
        dat <- rbind(dat[1:12,],dat[27:40,])
    }
    dat$Pos <- 1:nrow(dat)
    return(dat)
}

readTestSet29FromDamage_lik <- function(forward=TRUE){
    #This dataset was created in following fashion
    #nice -n 18 ./damage_lik -empirical -counts -includeSelf -i cgg10029_CTTGTA_L003_R1_.trunc_Eq_nucl_Onlynoseed_Mkdup_noCloseHuman_10K.sam -ref /home/hakon/data/4Hakon/Equus_caballus_genome_chrid.fasta 
    dat <- read.table("horse_dataset/dataset_29_from_dam_lik.txt",header=TRUE)
    if (forward){
        dat <- dat[1:12,]
    }else {
        dat <- rbind(dat[1:12,],dat[29:40,])
    }
    dat$Pos <- 1:nrow(dat)
    rownames(dat) <- 1:nrow(dat)
    return(dat)
}

readPhilSim <- function(forward=TRUE){
    #See /home/hakon/LocalWork/paraDamage/theModel/Phil_sim/README for simulation  
    dat <- read.table("/home/hakon/LocalWork/paraDamage/theModel/Phil_sim/data.txt",header=TRUE)
    dat[,"A"] <- dat[,"A"]+dat[,"A.C"]+dat[,"A.G"]+dat[,"A.T"]
    dat[,"C"] <- dat[,"C"]+dat[,"C.A"]+dat[,"C.G"]+dat[,"C.T"]
    dat[,"G"] <- dat[,"G"]+dat[,"G.A"]+dat[,"G.C"]+dat[,"G.T"]
    dat[,"T"] <- dat[,"T"]+dat[,"T.A"]+dat[,"T.C"]+dat[,"T.G"]
    if (forward){
        dat <- dat[1:12,]
    }else {
        dat <- rbind(dat[1:20,],dat[21:40,])
    }
    dat$Pos <- 1:nrow(dat)
    rownames(dat) <- 1:nrow(dat)
    return(dat)
}
