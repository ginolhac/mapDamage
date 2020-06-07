import csv

import mapdamage.seqtk as seqtk


def write_base_comp(fasta, destination):
    """Calculates the total base composition across all sequences in 'fasta'
    and writes them to 'destination' as CSV.
    """
    bases = {"A": 0, "C": 0, "G": 0, "T": 0}
    for stats in seqtk.comp(str(fasta)):
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
    """Read the base compition from a file created by write_base_comp
    """
    with open(filename, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            return row

    raise csv.Error("No rows found in %r" % (filename,))
