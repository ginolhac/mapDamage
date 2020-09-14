source("common.r")


OPT.LGDIST <- getArgument("LGDIST")
OPT.PDFOUT <- getArgument("PDFOUT")
OPT.MISINCORP <- getArgument("MISINCORP")
OPT.TITLE <- getArgument("TITLE")
OPT.VERSION <- getArgument("VERSION")

MISMATCHES <- c("C>T", "G>A")


plot.length.distribution <- function(tbl) {
  tbl <- aggregate(tbl$Occurences, tbl[, c("Kind", "Std", "Length")], sum)

  row <- 1
  data <- matrix(0, nrow = 4, ncol = max(tbl$Length))
  for (kind in c("se", "pe")) {
    for (strand in c("+", "-")) {
      # PE reads with unknown template lengths (Length = 0) are counted but not plotted
      subtbl <- subset(tbl, Std == strand & Kind == kind & Length > 0)

      data[row, subtbl$Length] <- subtbl$x
      row <- row + 1
    }
  }

  # The sum of space and width defines the width of each bar + padding
  space <- 0.2
  width <- 1.0

  # Start and end the axis on a number divisible by 10
  min_len <- floor(min(tbl$Length) / 10) * 10
  max_len <- ceiling(max(tbl$Length) / 10) * 10

  colors <- c(rgb(1:0, 0, 0:1, 1 / 2), grey.colors(3))
  barplot(data, width = width, space = space, border = NA, col = colors, axes = FALSE, axisnames = FALSE, main = "Length distribution", xlim = c(min_len, max_len) * (space + width))

  legend("topright",
    c("+ strand (SE)", "- strand (SE)", "+ strand (PE)", "- strand (PE)"),
    fill = colors, bty = "n", border = NA
  )

  mtext("Occurences", side = 2, line = 2.5, cex = 0.7)
  mtext("Read length", side = 1, line = 2, cex = 0.7)
  xcoord <- seq(min_len, max_len, 10)
  axis(side = 1, labels = xcoord, at = (xcoord - 0.5) * (width + space), las = 2, cex.axis = 0.6)
  axis(side = 2, labels = TRUE, las = 2, cex.axis = 0.6)
}


plot.cumulative.mutations <- function(tbl, end, mut, sid) {
  tbl <- aggregate(tbl[, MISMATCHES], tbl[, c("End", "Std", "Pos")], sum)
  subplus <- subset(tbl, Std == "+" & End == end)
  subminus <- subset(tbl, Std == "-" & End == end)
  plot(c(0, cumsum(subplus[, mut] / sum(subplus[, mut]))),
    type = "l", col = rgb(1, 0, 0, 1 / 2), lwd = 2, axes = FALSE
  )
  lines(c(0, cumsum(subminus[, mut] / sum(subminus[, mut]))),
    col = rgb(0, 0, 1, 1 / 2), lwd = 2
  )
  axis(side = 1, labels = TRUE, las = 2, cex.axis = 0.6)
  axis(side = sid, labels = seq(0, 1, 0.1), at = seq(0, 1, 0.1), las = 2, cex.axis = 0.6)
  mtext(mut, side = 3, line = 2, cex = 0.8)
  mtext("Read position", side = 1, line = 1.8, cex = 0.7)
  mtext("Cumulative frequencies", side = sid, line = 2.5, cex = 0.7)
  legend("topleft", c("+ strand", "- strand"),
    fill = rgb(1:0, 0, 0:1, 0.4), bty = "n",
    border = NA
  )
}


plot.library <- function(lengths, mut, title = "", subtitle = "") {
  par(oma = c(4, 2, 2, 2), mar = c(1, 2, 1, 2))
  layout(matrix(c(
    # Title
    1, 1,
    # lengths
    2, 2,
    # Cumulative mutation
    3, 4
  ),
  3, 2,
  byrow = TRUE
  ),
  heights = c(3, 20, 20)
  )

  # Plot title
  plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  mtext(title, 3, 0.5, cex = 1.3)
  mtext(subtitle, 3, -1.0, cex = 0.75)

  # Base compositions
  plot.length.distribution(lengths)

  par(mar = c(1, 2, 7, 1))
  plot.cumulative.mutations(mut, "5p", "C>T", 2)
  par(mar = c(1, 1, 7, 2))
  plot.cumulative.mutations(mut, "3p", "G>A", 4)
}


main <- function() {
  lengths <- read.mapDamage.table(OPT.LGDIST)
  misincorp <- read.mapDamage.table(OPT.MISINCORP)

  if (nrow(lengths) == 0) {
    write("No length distributions are available; cannot plot lengths!", stderr())
    return(0)
  }

  pdf(file = OPT.PDFOUT, title = paste("mapDamage-", OPT.VERSION, sep = ""))

  # Summary plot for entire BAM file
  plot.library(lengths, misincorp, title = OPT.TITLE)

  # Per-libraries plots (only if there are more than 1 library)
  for (lib in iterLibraries(lengths, min_libraries = 2)) {
    lib_len <- subset(lengths, Sample == lib$Sample & Library == lib$Library)
    lib_mis <- subset(misincorp, Sample == lib$Sample & Library == lib$Library)

    subtitle <- sprintf("Sample: %s, Library: %s", lib$Sample, lib$Library)
    plot.library(lib_len, lib_mis, title = OPT.TITLE, subtitle = subtitle)
  }

  # graphics.off() calls dev.off(), but doesn't warn about the null device
  graphics.off()

  return(0)
}

quit(status = main())
