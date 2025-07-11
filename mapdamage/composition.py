import csv

import mapdamage
import mapdamage.seqtk


def count_ref_comp(read, chrom, before, after, comp):
    """record basae composition in external genomic regions"""
    std = "-" if read.is_reverse else "+"

    _update_table(comp[chrom]["5p"][std], before, range(-len(before), 0))
    _update_table(comp[chrom]["3p"][std], after, range(1, len(after) + 1))


def count_read_comp(read, chrom, length, comp):
    """record base composition of read, discard marked nucleotides"""
    std, seq = "+", read.query
    if read.is_reverse:
        std, seq = "-", mapdamage.seq.revcomp(seq)

    _update_table(comp[chrom]["5p"][std], seq, range(1, length + 1))
    _update_table(comp[chrom]["3p"][std], reversed(seq), range(-1, -length - 1, -1))


def _update_table(table, sequence, indices):
    for index, nt in zip(indices, sequence):
        if nt in table:
            table[nt][index] += 1


def write_base_comp(fasta, destination):
    """Calculates the total base composition across all sequences in 'fasta'
    and writes them to 'destination' as CSV.
    """
    bases = {"A": 0, "C": 0, "G": 0, "T": 0}
    for stats in mapdamage.seqtk.comp(str(fasta)):
        for key in bases:
            bases[key] += stats[key]

    # calculate the base frequencies
    ba_su = sum(bases.values())
    for key in bases:
        bases[key] = bases[key] / ba_su

    with open(destination, "wt", newline="") as handle:
        writer = csv.writer(handle)

        header = ["A", "C", "G", "T"]
        writer.writerow(header)
        writer.writerow(bases[key] for key in header)


def read_base_comp(filename):
    """Read the base compition from a file created by write_base_comp"""
    with open(filename, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            return row

    raise csv.Error("No rows found in %r" % (filename,))
