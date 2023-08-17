source("common.r")

NUCLEOTIDES <- c("A", "C", "G", "T", "Total")
MISMATCHES <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
INSERTIONS <- c("->A", "->C", "->G", "->T")
DELETIONS <- c("A>-", "C>-", "G>-", "T>-")
CLIPPING <- c("S")
EVERYTHING <- c(NUCLEOTIDES, MISMATCHES, INSERTIONS, DELETIONS, CLIPPING)

OPT.COMP <- getArgument("COMP")
OPT.PDFOUT <- getArgument("PDFOUT")
OPT.AROUND <- getArgument("AROUND", as.integer)
OPT.MISINCORP <- getArgument("MISINCORP")
OPT.LENGTH <- getArgument("LENGTH", as.integer)
OPT.YMAX <- getArgument("YMAX", as.numeric)
OPT.FOLDER <- getArgument("FOLDER")
OPT.TITLE <- getArgument("TITLE")
OPT.VERSION <- getArgument("VERSION")


draw.open.rect <- function(xleft, ybottom, xright, ytop, padding = 0) {
  if (xleft < 0) {
    xpoints <- c(xleft, xright + padding, xright + padding, xleft)
  } else {
    xpoints <- c(xright, xleft - padding, xleft - padding, xright)
  }

  lines(xpoints, c(ytop, ytop, ybottom, ybottom), col = "darkgrey")
}


plot.base.composition <- function(tbl, base, color, around, ylabels.at = c(), xlabels = FALSE) {
  xcoords <- c(-around:-1, 1:around)

  plot.axis <- function(yaxis.at) {
    axis(side = yaxis.at, labels = (yaxis.at == ylabels.at), line = 0, las = 2, cex.axis = 0.8)
    if ((yaxis.at == 2) && any(yaxis.at == ylabels.at)) {
      mtext("Frequency", side = 2, line = 2.5, cex = 0.6)
    }

    if (xlabels) {
      axis(side = 1, labels = xcoords, at = xcoords, las = 2, cex.axis = 0.6)
    } else {
      axis(side = 1, labels = FALSE, at = xcoords)
    }
  }

  plot.frequencies <- function(end) {
    subtbl <- subset(tbl, End == end)
    plot(subtbl$Pos, subtbl[, base] / subtbl$Total,
      pch = ".",
      xlim = c(-around, around), ylim = c(0, 0.5),
      col = color, main = base, cex.axis = 0.8, las = 2,
      xlab = "", ylab = "", lab = c(2 * around, 6, 0.2),
      axes = FALSE
    )

    ycoords <- NULL
    for (i in xcoords) {
      ycoords <- append(ycoords, mean(subtbl[(subtbl$Pos == i), base] / subtbl$Total[(subtbl$Pos == i)], na.rm = T))
    }
    points(xcoords, ycoords, pch = 20, col = color, type = "b")
  }

  # 5p end
  par(mar = c(1, 2, 1, 0.25))
  plot.frequencies("5p")
  plot.axis(2)
  draw.open.rect(0, 0, around, 0.5)

  # 3p end
  par(mar = c(1, 0.25, 1, 2))
  plot.frequencies("3p")
  plot.axis(4)
  draw.open.rect(-around, 0, 0, 0.5)

  par(mar = c(1, 2, 1, 2))
}


calculate.mutation.table <- function(tbl, length) {
  tbl <- aggregate(tbl[, EVERYTHING], tbl[, c("End", "Pos")], sum)
  for (mismatch in MISMATCHES) {
    tbl[, mismatch] <- tbl[, mismatch] / tbl[, substr(mismatch, 1, 1)]
  }

  for (mismatch in c(INSERTIONS, DELETIONS, CLIPPING)) {
    tbl[, mismatch] <- tbl[, mismatch] / tbl[, "Total"]
  }

  return(tbl[tbl$Pos <= length, ])
}


plot.mutations <- function(mut, end, axis.at, start.i, end.i, modifier) {
  do.plot <- function(tbl, end, modifier, mismatches, color, width) {
    subtable <- tbl[tbl$End == end, c("Pos", mismatches)]
    rates <- rowSums(subtable) - subtable$Pos
    subtable <- aggregate(list(Rate = rates), list(Pos = subtable$Pos), sum)

    lines(subtable$Pos * modifier, subtable$Rate,
      xlim = c(1, OPT.LENGTH), ylim = c(0, OPT.YMAX),
      col = color, lwd = width
    )
  }

  plot(NA, xlim = c(start.i, end.i), ylim = c(0, OPT.YMAX), col = "grey", lwd = 1, type = "l", xlab = "", ylab = "", axes = FALSE)
  axis(side = 1, labels = start.i:end.i, at = start.i:end.i, las = 2, cex.axis = 0.8)
  axis(side = axis.at, labels = TRUE, las = 2, cex.axis = 0.8)
  if (end == "5p") {
    mtext("Frequency", side = 2, line = 2.5, cex = 0.6)
  }
  draw.open.rect(start.i, OPT.YMAX, end.i, -0.01, padding = 0.5)

  for (mismatch in MISMATCHES) {
    do.plot(mut, end, modifier, mismatch, "grey", 1)
  }

  do.plot(mut, end, modifier, CLIPPING, "orange", 1)
  do.plot(mut, end, modifier, DELETIONS, "green", 1)
  do.plot(mut, end, modifier, INSERTIONS, "purple", 1)
  do.plot(mut, end, modifier, "G>A", "blue", 2)
  do.plot(mut, end, modifier, "C>T", "red", 2)
}


plot.main <- function(mut, com, title = "", subtitle = "") {
  mut <- calculate.mutation.table(mut, OPT.LENGTH)
  com <- aggregate(com[, NUCLEOTIDES], com[, c("End", "Pos")], sum)

  par(oma = c(4, 2, 2, 2), mar = c(1, 2, 1, 2))
  layout(matrix(c(
    1, 1, 1, 1, # Title
    2, 3, 4, 5, # A,  C
    6, 7, 8, 9, # G,  T
    10, 10, 11, 11
  ), # Mismatches 5p, 3p
  4, 4,
  byrow = TRUE
  ),
  heights = c(3, 20, 20, 20)
  )

  # Plot title
  plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  mtext(title, 3, 0.5, cex = 1.3)
  mtext(subtitle, 3, -1.0, cex = 0.75)

  # Base compositions
  plot.base.composition(com, "A", "blue", OPT.AROUND, xlabels = FALSE, ylabels.at = c(2, 4))
  plot.base.composition(com, "C", "green", OPT.AROUND, xlabels = FALSE, ylabels.at = c(4))
  plot.base.composition(com, "G", "black", OPT.AROUND, xlabels = TRUE, ylabels.at = c(2, 4))
  plot.base.composition(com, "T", "red", OPT.AROUND, xlabels = TRUE, ylabels.at = c(4))

  # Misincorporation patterns
  par(mar = c(1, 2, 1, 1))
  plot.mutations(mut, "5p", 2, 1, OPT.LENGTH, 1)

  par(mar = c(1, 1, 1, 2))
  plot.mutations(mut, "3p", 4, -OPT.LENGTH, -1, -1)
}


main <- function() {
  pdf(file = OPT.PDFOUT, title = paste("mapDamage-", OPT.VERSION, " plot"))

  mut <- read.mapDamage.table(OPT.MISINCORP)
  com <- read.mapDamage.table(OPT.COMP)

  # Summary plot for entire BAM file
  plot.main(mut, com, title = OPT.TITLE)

  # Per-libraries plots (only if there are more than 1 library)
  for (lib in iterLibraries(mut, min_libraries = 2)) {
    lib_mut <- subset(mut, Sample == lib$Sample & Library == lib$Library)
    lib_com <- subset(com, Sample == lib$Sample & Library == lib$Library)

    subtitle <- sprintf("Sample: %s, Library: %s", lib$Sample, lib$Library)
    plot.main(lib_mut, lib_com, title = OPT.TITLE, subtitle = subtitle)
  }

  # graphics.off() calls dev.off(), but doesn't warn about the null device
  graphics.off()

  return(0)
}

quit(status = main())
