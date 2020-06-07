#Functions to load the raw data
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
        raw_dat_plus <- raw_dat[raw_dat$End=="5p" & raw_dat$Std=="+",] 
        raw_dat_minus <- raw_dat[raw_dat$End=="5p" & raw_dat$Std=="-",]
    } else {
        raw_dat_plus <- raw_dat[raw_dat$End=="3p" & raw_dat$Std=="+",]
        raw_dat_minus <- raw_dat[raw_dat$End=="3p" & raw_dat$Std=="-",] 
    }
    dat <- matrix(nrow=max(raw_dat_plus$Pos),ncol=length(colnames(raw_dat_plus))-3)
    colnames(dat) <- colnames(raw_dat_plus)[c(-1,-2,-3)]
    dat[,"Pos"] <- seq(from=1,to=nrow(dat),by=1)
    if (forward!=1){
        dat[,"Pos"] <- - dat[,"Pos"]
    }
    for (i in colnames(raw_dat)[c(-1,-2,-3,-4)]){
        dat[,i] <- sumTheChr(raw_dat_plus,i)+sumTheChr(raw_dat_minus,i)
    }
    return(dat)
}

joinFowAndRev <- function(fo,re,nrPos){
    #Joins the 5' and 3' end with nrPos bases for each end
    te <- fo[1:nrPos,]
    te2 <- re[nrPos:1,]
    out <- rbind(te,te2)
    out[,"Pos"] <- c(1:nrPos,-c(nrPos:1))
    return(out)
}

readBaseFreqs <- function(fol){
    #Get the nucleotide composition of the genome
    fil <- "dnacomp_genome.csv"
    dat <- read.csv(paste(fol,fil,sep=""),header=TRUE)
    return(c(dat$A,dat$C,dat$G,dat$T))
}
