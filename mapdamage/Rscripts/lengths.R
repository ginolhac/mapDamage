# Enable full backtraces on errors
on_error <- function(e)
  {
    traceback(2)
    quit(status = 1)
  }
options(error = on_error)


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


plot.lengthStd <- function(tbl, title) {
    data <- NULL
    min_len <- min(tbl$Length)
    max_len <- max(tbl$Length)
    for (kind in c("se", "pe")) {
        for (strand in c("+", "-")) {
            row <- NULL

            subtbl <- subset(tbl, Std == strand & Kind == kind)
            for (length in min_len:max_len) {
                row <- c(row, sum(subtbl[subtbl$Length == length, "Occurences"]))
            }

            data <- rbind(data, row)
        }
    }

    colnames(data) <- min_len:max_len

    colors <- c(rgb(1:0, 0, 0:1, 1/2), grey.colors(3))
    barplot(data, border = NA, col = colors, axes = FALSE, axisnames = FALSE, main = title)

    legend('topright',
           c('+ strand (SE)','- strand (SE)', '+ strand (PE)','- strand (PE)'),
           fill = colors, bty = 'n', border = NA)

    mtext("Occurences", side = 2, line = 2.5, cex = 0.7)
    mtext("Read length", side = 1, line = 2, cex = 0.7)
    xcoord = seq(min_len, max_len, 10)
    axis(side = 1, labels = xcoord, at = xcoord, las = 2, cex.axis = 0.6)
    axis(side = 2, labels = TRUE, las = 2, cex.axis = 0.6)
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

lg <- read.table(file = OPT.LGDIST, sep = "\t", header = TRUE, as.is = TRUE)
if(nrow(lg) == 0){
    write("No length distributions are available, plotting length distribution only works for single-end reads", stderr())
}else{

    pdf(file = OPT.PDFOUT, title = paste("mapDamage-", OPT.VERSION, sep=""))
    par(oma = c(4,2,2,2), mar = c(1,2,1,2))
    layout(matrix(c(1,1,  # Title
                    2,2,  # lengths
                    3,4), # Cumulative mutation
                    3, 2, byrow = TRUE),
           heights = c(3, 20, 20))
    
    # Plot title
    plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    mtext(OPT.TITLE, 3, cex = 1.3)
    
    # Base compositions
    plot.lengthStd(lg, "Length distribution")
    
    # Misincorporation patterns
    mut <- calculate.mutation.table(OPT.MISINCORP)
    
    par(mar = c(1, 2, 7, 1))
    plot.cumul.mutation(mut, "5p", "C>T", 2)
    par(mar = c(1, 1, 7, 2))
    plot.cumul.mutation(mut, "3p", "G>A", 4)
    
    # graphics.off() calls dev.off() for all devices but doesn't return anything (avoid null device message)
    graphics.off()
}
