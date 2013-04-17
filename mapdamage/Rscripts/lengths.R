args <- commandArgs(trailingOnly = TRUE)
OPT.LGDIST    <- args[1]
OPT.PDFOUT    <- args[2]
OPT.MISINCORP <- args[3]
OPT.LENGTH    <- args[4]
OPT.TITLE     <- args[5]
OPT.VERSION   <- args[6]

MISMATCHES  <- c("C>T", "G>A")

calculate.mutation.table <- function(filename){

    tbl <- read.table(file = filename, sep = "\t", header = TRUE, check.names = FALSE)
    tbl <- aggregate(tbl[, MISMATCHES], tbl[, c("End", "Std", "Pos")], sum)
    return(tbl)
}


plot.length <- function(tbl, title, color) {
    table <- aggregate(tbl$Occurences, by=list(tbl$Length), FUN=sum)
    names(table)<- c("Length", "Occurences")
    plot(table$Length, table$Occurences, type="h",
       col = color, main = title, cex.axis = 0.8, las = 2,
       xlab = "", ylab = "", axes = FALSE)
 
    mtext("Occurences", side = 2, line = 2.5, cex = 0.7)
    mtext("Read length", side = 1, line = 2, cex = 0.7)
    xcoord = seq(min(table$Length), max(table$Length), 10)
    axis(side = 1, labels = xcoord, at = xcoord, las = 2, cex.axis = 0.6)
    axis(side = 2, labels = TRUE, las = 2, cex.axis = 0.6)
}

plot.lengthStd <- function(tbl, title) {
    subplus = subset(tbl, Std == "+")
    subminus = subset(tbl, Std == "-")
    plot(subplus$Length, subplus$Occurences, type="h", col=rgb(1,0,0,1/2), main = title, axes = FALSE)
    lines(subminus$Length, subminus$Occurences, type="h", col=rgb(0,0,1,1/2))
    mtext("Occurences", side = 4, line = 2.5, cex = 0.7)
    mtext("Read length", side = 1, line = 2, cex = 0.7)
    xcoord = seq(min(subplus$Length), max(subplus$Length), 10)
    axis(side = 1, labels = xcoord, at = xcoord, las = 2, cex.axis = 0.6)
    axis(side = 4, labels = TRUE, las = 2, cex.axis = 0.6)
    legend('topright',c('+ strand','- strand'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
}

plot.cumul.mutation <- function(tbl, end, mut, sid) {
    subplus = subset(tbl, Std == "+" & End == end)
    subminus = subset(tbl, Std == "-" & End == end)
    plot(c(0, cumsum(subplus[, mut]/sum(subplus[, mut]))), 
	 type = "l", col = rgb(1,0,0,1/2), lwd = 2, axes = FALSE)
    lines(c(0, cumsum(subminus[, mut]/sum(subminus[, mut]))), 
	 col = rgb(0,0,1,1/2), lwd = 2)
    axis(side = 1, labels = TRUE, las = 2, cex.axis = 0.6)
    axis(side = sid, labels = seq(0,1,0.1), at = seq(0,1,0.1), las = 2, cex.axis = 0.6)
    mtext(mut, side = 3, line = 2, cex = 0.8)
    mtext("Read position", side = 1, line = 1.8, cex = 0.7)
    mtext("Cumulative frequencies", side = sid, line = 2.5, cex = 0.7)
    legend('topleft',c('+ strand','- strand'),
       	  fill = rgb(1:0,0,0:1,0.4), bty = 'n',
          border = NA)
}


pdf(file = OPT.PDFOUT, title = paste("mapDamage-", OPT.VERSION, " plot"))
par(oma = c(4,2,2,2), mar = c(1,2,1,2))
layout(matrix(c(1,1,  # Title
                2,3,  # lengths
                4,5), # Cumulative mutation
                3, 2, byrow = TRUE),
       heights = c(3, 20, 20))

# Plot title
plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
mtext(OPT.TITLE, 3, cex = 1.3)

# Base compositions
lg <- read.table(file = OPT.LGDIST, sep = "\t", header = TRUE, as.is = TRUE)
plot.length(lg, "Single-end read length distribution", "black")
plot.lengthStd(lg, "Single-end read length per strand")

# Misincorporation patterns
mut <- calculate.mutation.table(OPT.MISINCORP)

par(mar = c(1, 2, 7, 1))
plot.cumul.mutation(mut, "5p", "C>T", 2)
par(mar = c(1, 1, 7, 2))
plot.cumul.mutation(mut, "3p", "G>A", 4)

# graphics.off() calls dev.off() for all devices but doesn't return anything (avoid null device message)
graphics.off()

