
readMapDamData <- function(folder, termini="both", sub_length=12) {
    nucleotides <- c("A", "C", "G", "T")
    mismatches <- c("A.C", "A.G", "A.T",
                    "C.A", "C.G", "C.T",
                    "G.A", "G.C", "G.T",
                    "T.A", "T.C", "T.G")

    data <- read.table(paste(folder, "misincorporation.txt", sep=""), header=TRUE)
    data <- data[data$Pos <= sub_length, ]
    data[data$End == "3p", "Pos"] <- -data[data$End == "3p", "Pos"]

    if (termini == "5p" || termini == "3p") {
        data <- data[data$End == termini,]
    } else if (termini != "both") {
        abort("invalid termini '%s'", termini)
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
