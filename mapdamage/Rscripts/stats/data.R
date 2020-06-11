
readMapDamData <- function(folder, direction="both", sub_length=12) {
    nucleotides <- c("A", "C", "G", "T")
    mismatches <- c("A.C", "A.G", "A.T",
                    "C.A", "C.G", "C.T",
                    "G.A", "G.C", "G.T",
                    "T.A", "T.C", "T.G")

    data <- read.table(paste(folder, "misincorporation.txt", sep=""), header=TRUE)
    data <- data[data$Pos <= sub_length, ]
    data[data$End == "3p", "Pos"] <- -data[data$End == "3p", "Pos"]

    if (direction == "forward") {
        data <- data[data$End == "5p",]
    } else if (direction == "reverse") {
        data <- data[data$End == "3p",]
    } else if (direction != "both") {
        abort("invalid direction '%s'", direction)
    }

    return(aggregate(
        data[, c(nucleotides, mismatches)],
        list(Pos=data[, "Pos"]),
        sum
    ))
}


readBaseFreqs <- function(folder) {
    # Get the nucleotide composition of the genome
    data <- read.csv(paste(folder, "dnacomp_genome.csv", sep=""), header=TRUE)

    return(c(data$A, data$C, data$G, data$T))
}
