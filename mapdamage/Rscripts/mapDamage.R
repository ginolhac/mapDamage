# R script mapdamage

args <- commandArgs(trailingOnly = TRUE)

#print(args)

comp<-args[1]
pdfout<-args[2]
around<-as.numeric(args[3])
misincorp<-args[4]
lg<-as.numeric(args[5])
ymax<-as.numeric(args[6])
folder<-args[7]
title<-args[8]
version<-args[9]

com<-read.table(file=comp,sep="\t", header=TRUE, as.is=TRUE)
five=subset(com, End=="5p")
three=subset(com, End=="3p")
pdf(file=pdfout, title=paste("mapDamage-",version," plot"))
par(oma=c(4,2,2,2),mar=c(1,2,1,1))
layout(matrix(c(1,2,3,4,5,6,7,8,9,9,10,10), 3, 4, byrow=TRUE))
# for letter A, before - read
n<-numeric
plot(five$Pos,five$A/five$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="blue",main="A",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=2,labels=TRUE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=FALSE)
mtext("Frequency",side=2,line=2.5,cex=0.7)
rect(0.5,0,around+0.5,0.5,border="darkgrey")
segments(around+0.5,0,around+0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(five$A[(five$Pos==i)]/five$Tot[(five$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="blue",type="b")
# for letter A, read - after
n<-numeric
plot(three$Pos,three$A/three$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="blue",main="A",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=4,labels=FALSE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=FALSE)
rect(-around-0.5,0,-0.5,0.5,border="darkgrey")
segments(-around-0.5,0,-around-0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(three$A[(three$Pos==i)]/three$Tot[(three$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="blue",type="b")							 
mtext(title, side=3, line=1.2, cex=0.8)

# for letter C, before - read
n<-numeric
plot(five$Pos,five$C/five$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="green",main="C",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=2,labels=FALSE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=FALSE)
rect(0.5,0,around+0.5,0.5,border="darkgrey")
segments(around+0.5,0,around+0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(five$C[(five$Pos==i)]/five$Tot[(five$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="green",type="b")
# for letter C, read - after
n<-numeric
plot(three$Pos,three$C/three$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="green",main="C",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=4,labels=TRUE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=FALSE)
rect(-around-0.5,0,-0.5,0.5,border="darkgrey")
segments(-around-0.5,0,-around-0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(three$C[(three$Pos==i)]/three$Tot[(three$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="green",type="b")							 

# for letter G, before - read
n<-numeric
plot(five$Pos,five$G/five$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="black",main="G",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=2,labels=TRUE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=TRUE, las=2, cex.axis=0.6)
mtext("Frequency",side=2,line=2.5,cex=0.7)
rect(0.5,0,around+0.5,0.5,border="darkgrey")
segments(around+0.5,0,around+0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(five$G[(five$Pos==i)]/five$Tot[(five$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="black",type="b")
# for letter G, read - after
n<-numeric
plot(three$Pos,three$G/three$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="black",main="G",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=4,labels=FALSE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=TRUE, las=2, cex.axis=0.6)
rect(-around-0.5,0,-0.5,0.5,border="darkgrey")
segments(-around-0.5,0,-around-0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(three$G[(three$Pos==i)]/three$Tot[(three$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="black",type="b")

# for letter T, before - read
n<-numeric
plot(five$Pos,five$T/five$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="red",main="T",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=2,labels=FALSE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=TRUE, las=2, cex.axis=0.6)
rect(0.5,0,around+0.5,0.5,border="darkgrey")
segments(around+0.5,0,around+0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(five$T[(five$Pos==i)]/five$Tot[(five$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="red",type="b")
# for letter T, read - after
n<-numeric
plot(three$Pos,three$T/three$Tot,pch='.',xlim=c(-around,around),ylim=c(0,0.5),col="red",main="T",cex.axis=0.8,las=2,xlab="",ylab="",lab=c(2*around,6,0.2), axes=FALSE)
axis(side=4,labels=TRUE,line=0,las=2,cex.axis=0.8)
axis(side=1,labels=TRUE, las=2, cex.axis=0.6)
rect(-around-0.5,0,-0.5,0.5,border="darkgrey")
segments(-around-0.5,0,-around-0.5,0.5,col="white",lwd=2)
for (i in c(-around:-1, 1:around)) { 
  n<-c(n,mean(three$T[(three$Pos==i)]/three$Tot[(three$Pos==i)],na.rm=T)) 
}
points(c(-around:-1,1:around),n[2:((2*around)+1)],pch=20,col="red",type="b")

# Misincorporation patterns
# first from 5'-ends
mut<-read.table(file=misincorp,sep="\t",header=TRUE, check.names = FALSE)
nucl<-subset(mut, End=="5p")
vec<-(1:lg)
sumhit<-data.frame(cbind(vec,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
insertion<-0
for (i in vec)	{
  sumhit[i,2]<-sum(nucl[(nucl$Pos == i), "G>A"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,3]<-sum(nucl[(nucl$Pos == i), "C>T"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,4]<-sum(nucl[(nucl$Pos == i), "A>G"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,5]<-sum(nucl[(nucl$Pos == i), "T>C"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,6]<-sum(nucl[(nucl$Pos == i), "A>C"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,7]<-sum(nucl[(nucl$Pos == i), "A>T"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,8]<-sum(nucl[(nucl$Pos == i), "C>G"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,9]<-sum(nucl[(nucl$Pos == i), "C>A"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,10]<-sum(nucl[(nucl$Pos == i), "T>G"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,11]<-sum(nucl[(nucl$Pos == i), "T>A"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,12]<-sum(nucl[(nucl$Pos == i), "G>C"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,13]<-sum(nucl[(nucl$Pos == i), "G>T"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,14]<-sum(nucl[(nucl$Pos == i), "A>-"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,15]<-sum(nucl[(nucl$Pos == i), "T>-"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,16]<-sum(nucl[(nucl$Pos == i), "C>-"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,17]<-sum(nucl[(nucl$Pos == i), "G>-"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,22]<-sum(nucl[(nucl$Pos == i), "S"])/sum(nucl[(nucl$Pos == i), "Tot"])
  NTs <- c("A", "C", "G", "T")
  for (j in seq(NTs)) {
    insertion[j]<-sprintf("->%s", NTs[j])
    sumhit[i,17 + j] <-sum(nucl[(nucl$Pos == i), insertion[j]]) /sum(nucl[(nucl$Pos == i), NTs])
  }
}
# write a table for C>T frequencies
write.table(sumhit[,3], file=paste(folder,"/5pCtoT_freq.txt",sep=""), sep="\t", quote=FALSE, col.names=c("pos\t5pC>T"))
plot(sumhit$vec,sumhit$V4,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1,type="l",xlab="",ylab="",cex.lab=0.8,cex.axis=0.6,las=2,lab=c(lg,11,0.1),axes=FALSE)
axis(side=1,labels=TRUE,las=2,cex.axis=0.8)
axis(side=2,labels=TRUE,las=2,cex.axis=0.8)
rect(0.5,-0.01,lg+0.5,ymax,border="darkgrey")
segments(lg+0.5,-0.01,lg+0.5,ymax,col="white",lwd=2)
mtext("Frequency",side=2,line=2.8,cex=0.7)
lines(sumhit$vec,sumhit$V5,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V6,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V7,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V8,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V9,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V10,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V11,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V12,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
lines(sumhit$vec,sumhit$V13,xlim=c(1,lg),ylim=c(0,ymax),col="grey",lwd=1)
# sofclipping in orange
lines(sumhit$vec,sumhit$V22,xlim=c(1,lg),ylim=c(0,ymax),col="orange",lwd=1)
# deletio in green
lines(sumhit$vec,(sumhit$V14+sumhit$V15+sumhit$V16+sumhit$V17),xlim=c(1,lg),ylim=c(0,ymax),col="green",lwd=1)
# insertions in pruple
lines(sumhit$vec,(sumhit$V18+sumhit$V19+sumhit$V20+sumhit$V21),xlim=c(1,lg),ylim=c(0,ymax),col="purple",lwd=1)
# G>A in blue
lines(sumhit$vec,sumhit$V2,xlim=c(1,lg),ylim=c(0,ymax),col="blue",lwd=2)
# C>T in red
lines(sumhit$vec,sumhit$V3,xlim=c(1,lg),ylim=c(0,ymax),col="red",lwd=2)

# then from 3'-ends
nucl=subset(mut, End=="3p")
vec<-(1:lg)
sumhit<-data.frame(cbind(vec,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
insertion<-0
for (i in vec)	{
  sumhit[i,2]<-sum(nucl[(nucl$Pos == i), "G>A"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,3]<-sum(nucl[(nucl$Pos == i), "C>T"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,4]<-sum(nucl[(nucl$Pos == i), "A>G"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,5]<-sum(nucl[(nucl$Pos == i), "T>C"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,6]<-sum(nucl[(nucl$Pos == i), "A>C"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,7]<-sum(nucl[(nucl$Pos == i), "A>T"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,8]<-sum(nucl[(nucl$Pos == i), "C>G"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,9]<-sum(nucl[(nucl$Pos == i), "C>A"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,10]<-sum(nucl[(nucl$Pos == i), "T>G"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,11]<-sum(nucl[(nucl$Pos == i), "T>A"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,12]<-sum(nucl[(nucl$Pos == i), "G>C"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,13]<-sum(nucl[(nucl$Pos == i), "G>T"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,14]<-sum(nucl[(nucl$Pos == i), "A>-"])/sum(nucl[(nucl$Pos == i), "A"])
  sumhit[i,15]<-sum(nucl[(nucl$Pos == i), "T>-"])/sum(nucl[(nucl$Pos == i), "T"])
  sumhit[i,16]<-sum(nucl[(nucl$Pos == i), "C>-"])/sum(nucl[(nucl$Pos == i), "C"])
  sumhit[i,17]<-sum(nucl[(nucl$Pos == i), "G>-"])/sum(nucl[(nucl$Pos == i), "G"])
  sumhit[i,22]<-sum(nucl[(nucl$Pos == i), "S"])/sum(nucl[(nucl$Pos == i), "Tot"])
  NTs <- c("A", "C", "G", "T")
  for (j in seq(NTs)) {
    insertion[j]<-sprintf("->%s", NTs[j])
    sumhit[i,17 + j] <-sum(nucl[(nucl$Pos == i), insertion[j]]) /sum(nucl[(nucl$Pos == i), NTs])
  }
}
#  write a table for G>A frequencies
write.table(sumhit[,2], file=paste(folder,"/5pGtoA_freq.txt",sep=""), sep="\t", quote=FALSE, col.names=c("pos\t5pG>A"))
plot(-sumhit$vec,sumhit$V4,xlim=c(-lg, -1),ylim=c(0,ymax),col="grey",lwd=1,type="l",xlab="",ylab="",cex.lab=0.8,cex.axis=0.6,las=2,lab=c(lg,11,0.1),axes=FALSE)
axis(side=1,labels=TRUE,las=2,cex.axis=0.8)
axis(side=4,labels=TRUE,las=2,cex.axis=0.8)
rect(-lg-0.5,-0.01,-0.5,ymax,border="darkgrey")
segments(-lg-0.5,-0.01,-lg-0.5,ymax,col="white",lwd=2)
lines(-sumhit$vec,sumhit$V5,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V6,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V7,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V8,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V9,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V10,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V11,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V12,col="grey",lwd=1)
lines(-sumhit$vec,sumhit$V13,col="grey",lwd=1)
# sofclipping in orange
lines(-sumhit$vec,sumhit$V22,col="orange",lwd=1)
# deletions in green
lines(-sumhit$vec,(sumhit$V14+sumhit$V15+sumhit$V16+sumhit$V17),col="green",lwd=1)
# insertions in pruple
lines(-sumhit$vec,(sumhit$V18+sumhit$V19+sumhit$V20+sumhit$V21),col="purple",lwd=1)
# G>A in blue
lines(-sumhit$vec,sumhit$V2,col="blue",lwd=2)
# C>T in red
lines(-sumhit$vec,sumhit$V3,ylim=c(0,ymax),col="red",lwd=2)

# graphics.off() calls dev.off() for all devices but doesn't return anything (avoid null device message)
graphics.off()

