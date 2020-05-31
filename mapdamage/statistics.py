import collections
import logging
import os

import mapdamage


class MisincorporationRates:
    def __init__(self, length):
        self.length = length
        self.data = {}
        for end in ("5p", "3p"):
            table_end = self.data[end] = {}
            for strand in ("+", "-"):
                table_strand = table_end[strand] = {}
                for mut in mapdamage.seq.HEADER:
                    table_strand[mut] = dict.fromkeys(range(length), 0)

    def update(self, read, seq, refseq, end):
        strand = "-" if read.is_reverse else "+"
        subtable = self.data[end][strand]

        for index, nt_seq, nt_ref in zip(range(self.length), seq, refseq):
            if nt_seq in "ACGT-" and nt_ref in "ACGT-":
                if nt_ref != "-":
                    # record base composition in the reference, only A, C, G, T
                    subtable[nt_ref][index] += 1

                # Most ref/seq pairs will be identical
                if nt_ref != nt_seq:
                    mut = "%s>%s" % (nt_ref, nt_seq)
                    subtable[mut][index] += 1

    def update_soft_clipping(self, read):
        def update_table(end, std, bases):
            for i in range(0, min(bases, self.length)):
                self.data[end][std]["S"][i] += 1

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
        with open(filepath, "wt") as handle:
            _write_freq_table(self.data, mapdamage.seq.HEADER, handle, offset=1)


class DNAComposition:
    def __init__(self, around, length):
        keys = {
            "3p": list(range(-length, 0)) + list(range(1, around + 1)),
            "5p": list(range(-around, 0)) + list(range(1, length + 1)),
        }

        self.data = {}
        for end in ("5p", "3p"):
            tab_end = self.data[end] = {}
            for strand in ("+", "-"):
                tab_strand = tab_end[strand] = {}
                for letters in mapdamage.seq.LETTERS:
                    tab_strand[letters] = dict.fromkeys(keys[end], 0)

    def update_read(self, read, length):
        strand, seq = "+", read.query
        if read.is_reverse:
            strand, seq = "-", mapdamage.seq.revcomp(seq)

        self._update_table(self.data["5p"][strand], seq, range(1, length + 1))
        self._update_table(
            self.data["3p"][strand], reversed(seq), range(-1, -length - 1, -1)
        )

    def update_reference(self, read, before, after):
        strand = "-" if read.is_reverse else "+"

        self._update_table(self.data["5p"][strand], before, range(-len(before), 0))
        self._update_table(self.data["3p"][strand], after, range(1, len(after) + 1))

    def write(self, filepath):
        with open(filepath, "wt") as handle:
            columns = mapdamage.seq.LETTERS + ("Total",)
            _write_freq_table(self.data, columns, handle)

    def _update_table(self, table, sequence, indices):
        for index, nt in zip(indices, sequence):
            if nt in table:
                table[nt][index] += 1


class FragmentLengths:
    def __init__(self):
        self.data = {
            (kind, strand): collections.defaultdict(int)
            for kind in ("pe", "se")
            for strand in ("+", "-")
        }

    def update(self, read, coordinate):
        strand = "-" if read.is_reverse else "+"

        if read.is_paired:
            if read.is_read1 and read.is_proper_pair:
                length = abs(read.template_length)
                self.data[("pe", strand)][length] += 1
        else:
            length = max(coordinate) - min(coordinate)
            self.data[("se", strand)][length] += 1

    def write(self, filepath):
        with open(filepath, "wt") as handle:
            handle.write("Std\tKind\tLength\tOccurences\n")
            for (pe_or_se, strand), lengths in sorted(self.data.items()):
                for length, count in sorted(lengths.items()):
                    handle.write("%s\t%s\t%d\t%d\n" % (strand, pe_or_se, length, count))


def check_table_and_warn_if_dmg_freq_is_low(folder):
    """ Returns true if the damage frequencies are too low to allow
    Bayesian estimation of DNA damages, i.e < 1% at first position.
    """
    total = 0.0
    logger = logging.getLogger(__name__)
    for filename in ("5pCtoT_freq.txt", "3pGtoA_freq.txt"):
        if not os.path.exists(os.path.join(folder, filename)):
            logger.error(
                "Required table has not been created (%r), bayesian computation cannot be performed",
                filename,
            )
            return False

        with open(os.path.join(folder, filename)) as handle:
            for line in handle:
                freq = line.strip().split("\t")
                if freq[0] == "1":
                    total += float(freq[1])
                    break
            else:
                logger.error(
                    "Could not find pos = 1 in table %r, bayesian computation cannot be performed",
                    filename,
                )
                return False

    if total < 0.01:
        logger.warning(
            "DNA damage levels are too low, the Bayesian computation should not be performed (%f < 0.01)",
            total,
        )

    return True


def _write_freq_table(table, columns, out, offset=0):
    out.write("Chr\tEnd\tStd\tPos\t%s\n" % ("\t".join(columns)))

    for (end, strands) in sorted(table.items()):
        for (strand, subtable) in sorted(strands.items()):
            subtable["Total"] = {}
            for index in sorted(subtable[columns[0]]):
                subtable["Total"][index] = sum(
                    subtable[letter][index] for letter in mapdamage.seq.LETTERS
                )

                out.write("*\t%s\t%s\t%d" % (end, strand, index + offset))
                for base in columns:
                    out.write("\t%d" % subtable[base][index])
                out.write("\n")
