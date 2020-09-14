import collections
import csv
import logging
import os

import mapdamage


class MisincorporationRates:
    def __init__(self, libraries, length):
        self.length = length
        self.data = {}
        for library in libraries:
            table_lib = self.data[library] = {}
            for end in ("5p", "3p"):
                table_end = table_lib[end] = {}
                for strand in ("+", "-"):
                    table_strand = table_end[strand] = {}
                    for mut in mapdamage.seq.HEADER:
                        table_strand[mut] = dict.fromkeys(range(length), 0)

    def update(self, read, seq, refseq, end, library):
        strand = "-" if read.is_reverse else "+"
        subtable = self.data[library][end][strand]

        for index, nt_seq, nt_ref in zip(range(self.length), seq, refseq):
            if nt_seq in "ACGT-" and nt_ref in "ACGT-":
                if nt_ref != "-":
                    # record base composition in the reference, only A, C, G, T
                    subtable[nt_ref][index] += 1

                # Most ref/seq pairs will be identical
                if nt_ref != nt_seq:
                    mut = "%s>%s" % (nt_ref, nt_seq)
                    subtable[mut][index] += 1

    def update_soft_clipping(self, read, library):
        def update_table(end, std, bases):
            for i in range(0, min(bases, self.length)):
                self.data[library][end][std]["S"][i] += 1

        strand = "-" if read.is_reverse else "+"
        for nbases, idx in mapdamage.align.parse_cigar(read.cigar, 4):
            if idx == 0:
                # Soft-clipping at the left side of the alignment
                end = "3p" if read.is_reverse else "5p"
            else:
                # Soft-clipping at the right side of the alignment
                end = "5p" if read.is_reverse else "3p"

            update_table(end, strand, nbases)

    def write(self, filepath):
        with filepath.open("wt") as handle:
            _write_freq_table(self.data, mapdamage.seq.HEADER, handle, offset=1)


class DNAComposition:
    def __init__(self, libraries, around, length):
        keys = {
            "3p": list(range(-length, 0)) + list(range(1, around + 1)),
            "5p": list(range(-around, 0)) + list(range(1, length + 1)),
        }

        self.data = {}
        for library in libraries:
            tab_lib = self.data[library] = {}
            for end in ("5p", "3p"):
                tab_end = tab_lib[end] = {}
                for strand in ("+", "-"):
                    tab_strand = tab_end[strand] = {}
                    for letters in mapdamage.seq.LETTERS:
                        tab_strand[letters] = dict.fromkeys(keys[end], 0)

    def update_read(self, read, length, library):
        strand, seq = "+", read.query
        if read.is_reverse:
            strand, seq = "-", mapdamage.seq.revcomp(seq)

        self._update_table(self.data[library]["5p"][strand], seq, range(1, length + 1))
        self._update_table(
            self.data[library]["3p"][strand], reversed(seq), range(-1, -length - 1, -1)
        )

    def update_reference(self, read, before, after, library):
        strand = "-" if read.is_reverse else "+"

        self._update_table(
            self.data[library]["5p"][strand], before, range(-len(before), 0)
        )
        self._update_table(
            self.data[library]["3p"][strand], after, range(1, len(after) + 1)
        )

    def write(self, filepath):
        with filepath.open("wt") as handle:
            columns = mapdamage.seq.LETTERS + ("Total",)
            _write_freq_table(self.data, columns, handle)

    def _update_table(self, table, sequence, indices):
        for index, nt in zip(indices, sequence):
            if nt in table:
                table[nt][index] += 1


class FragmentLengths:
    def __init__(self, libraries):
        self.data = {
            library: {
                (kind, strand): collections.defaultdict(int)
                for kind in ("pe", "se")
                for strand in ("+", "-")
            }
            for library in libraries
        }

    def update(self, read, library):
        strand = "-" if read.is_reverse else "+"
        data = self.data[library]

        if read.is_paired:
            if read.is_read1 and read.is_proper_pair:
                length = abs(read.template_length)
                data[("pe", strand)][length] += 1
        else:
            data[("se", strand)][read.reference_length] += 1

    def write(self, filepath):
        with open(filepath, "wt") as handle:
            handle.write("Sample\tLibrary\tStd\tKind\tLength\tOccurences\n")
            for (sample, library), reads in sorted(self.data.items()):
                for (pe_or_se, strand), lengths in sorted(reads.items()):
                    for length, count in sorted(lengths.items()):
                        handle.write(
                            "%s\t%s\t%s\t%s\t%d\t%d\n"
                            % (sample, library, strand, pe_or_se, length, count)
                        )


def check_table_and_warn_if_dmg_freq_is_low(folder):
    """Returns true if the damage frequencies are too low to allow
    Bayesian estimation of DNA damages, i.e < 1% at first position.
    """
    logger = logging.getLogger(__name__)
    filename = "misincorporation.txt"
    mismatches = {
        "5p": {"C": 0, "C>T": 0},
        "3p": {"G": 0, "G>A": 0},
    }

    try:
        with open(os.path.join(folder, filename), newline="") as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            if not reader.fieldnames:
                logger.error("%r is empty; please re-run mapDamage", filename)
                return False

            for row in reader:
                if int(row["Pos"]) == 1:
                    counts = mismatches[row["End"]]
                    for key in counts:
                        counts[key] += int(row[key])
    except (csv.Error, IOError, OSError, KeyError) as error:
        logger.error("Error reading misincorporation table: %s", error)
        return False

    if not (mismatches["5p"]["C"] and mismatches["3p"]["G"]):
        logger.error(
            "Insufficient data in %r; cannot perform Bayesian computation", filename
        )
        return False

    total = 0.0
    total += mismatches["5p"]["C>T"] / mismatches["5p"]["C"]
    total += mismatches["3p"]["G>A"] / mismatches["3p"]["G"]

    if total < 0.01:
        logger.warning(
            "DNA damage levels are too low, the Bayesian computation should not be "
            "performed (%f < 0.01)",
            total,
        )

    return True


def _write_freq_table(table, columns, out, offset=0):
    out.write("Sample\tLibrary\tEnd\tStd\tPos\t%s\n" % ("\t".join(columns)))

    for (sample, library), ends in sorted(table.items()):
        for (end, strands) in sorted(ends.items()):
            for (strand, subtable) in sorted(strands.items()):
                subtable["Total"] = {}
                for index in sorted(subtable[columns[0]]):
                    subtable["Total"][index] = sum(
                        subtable[letter][index] for letter in mapdamage.seq.LETTERS
                    )

                    row = [sample, library, end, strand, str(index + offset)]
                    row.extend(str(subtable[base][index]) for base in columns)

                    out.write("\t".join(row))
                    out.write("\n")
