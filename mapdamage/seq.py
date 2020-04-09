import logging

# from Martin Kircher, to complement DNA
TABLE = str.maketrans("TGCAMRWSYKVHDBtgcamrwsykvhdb", "ACGTKYWSRMBDHVacgtkywsrmbdhv")

LETTERS = ("A", "C", "G", "T")
MUTATIONS = (
    "G>A",
    "C>T",
    "A>G",
    "T>C",
    "A>C",
    "A>T",
    "C>G",
    "C>A",
    "T>G",
    "T>A",
    "G>C",
    "G>T",
    "A>-",
    "T>-",
    "C>-",
    "G>-",
    "->A",
    "->T",
    "->C",
    "->G",
    "S",
)
HEADER = LETTERS + ("Total",) + MUTATIONS


def revcomp(seq):
    """ return reverse complemented string """
    return seq.translate(TABLE)[::-1]


def record_lg(read, coordinate, tab):
    """ record global length distribution
  don't record paired reads as they are normally not used for aDNA """
    std = "-" if read.is_reverse else "+"

    length = max(coordinate) - min(coordinate)
    if not read.is_paired:
        tab[std][length] = tab[std][length] + 1

    return tab


def read_fasta_index(filename):
    """ from a fasta index file, fai, return dictionary of references:lengths """
    logger = logging.getLogger(__name__)

    fai = {}
    with open(filename, "r") as handle:
        for line in handle:
            ref = line.split("\t")
            if len(ref) != 5:
                logger.error(
                    "Line %i in %r contains wrong number of fields, found %i, expected 5:",
                    line,
                    filename,
                    len(ref),
                )
                return None

            try:
                fai[ref[0]] = int(ref[1])
            except ValueError:
                logger.error(
                    "Length at line %i in %r is not a number; found %r",
                    line,
                    filename,
                    ref[1],
                )
                return None

    if not fai:
        logger.error("Error: Index for %r does contain any sequences.", filename)
        logger.error("Please ensure that FASTA file is valid, and")
        logger.error("re-index file using 'samtool faidx'.")
        return None

    return fai


def compare_sequence_dicts(fasta_dict, bam_dict):
    """Compares a FASTA and BAM sequence dictionary, and prints any differences.
  Returns true if all required sequences are found, false otherwise."""
    if fasta_dict == bam_dict:
        return True

    logger = logging.getLogger(__name__)
    common = set(fasta_dict) & set(bam_dict)
    if not common:
        logger.error("BAM and FASTA file have no sequence names in common")
        return False

    # Check that the lengths of required sequences match (fatal error)
    different = []
    for key in sorted(common):
        if fasta_dict[key] != bam_dict[key]:
            different.append((key, fasta_dict[key], bam_dict[key]))

    if different:
        logger.error("Length of required FASTA sequences differ:")
        for values in different:
            logger.error(" - %s: %i vs %i bp" % values)

    # Check for sequences only found in the BAM file (fatal errors)
    bam_only = set(bam_dict) - common
    if bam_only:
        logger.error("Sequences not found in FASTA:")
        for key in bam_only:
            logger.error("%s (%i bp)", key, bam_dict[key])

    # Check for sequences only found in the BAM file (fatal errors)
    fasta_only = set(fasta_dict) - common
    if fasta_only:
        logger.warning("FASTA file contains extra sequences:")
        for key in fasta_only:
            logger.warning(" - %s = %i bp" % (key, fasta_dict[key]))

    return not (different or bam_only)
