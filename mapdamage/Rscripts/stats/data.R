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
    if (forward!=1){
        dat[,"Pos"] <- - dat[,"Pos"]
    }
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
