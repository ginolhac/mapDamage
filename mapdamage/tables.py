import collections
import logging
import os

import mapdamage
from mapdamage.version import __version__


def initialize_mut(length):
    tab = {}
    for end in ("5p", "3p"):
        tab_end = tab[end] = {}
        for std in ("+", "-"):
            tab_std = tab_end[std] = {}
            for mut in mapdamage.seq.HEADER:
                tab_std[mut] = dict.fromkeys(range(length), 0)

    return tab


def print_mut(mut, opt, out):
    _print_freq_table(mut, mapdamage.seq.HEADER, opt, out, offset=1)


def initialize_comp(around, length):
    keys = {
        "3p": list(range(-length, 0)) + list(range(1, around + 1)),
        "5p": list(range(-around, 0)) + list(range(1, length + 1)),
    }

    tab = {}
    for end in ("5p", "3p"):
        tab_end = tab[end] = {}
        for std in ("+", "-"):
            tab_std = tab_end[std] = {}
            for letters in mapdamage.seq.LETTERS:
                tab_std[letters] = dict.fromkeys(keys[end], 0)

    return tab


def print_comp(comp, opt, out):
    columns = mapdamage.seq.LETTERS + ("Total",)
    _print_freq_table(comp, columns, opt, out)


def initialize_lg():
    return {
        (kind, strand): collections.defaultdict(int)
        for kind in ("pe", "se")
        for strand in ("+", "-")
    }


def print_lg(tab, opt, out):
    out.write("# table produced by mapDamage version %s\n" % __version__)
    out.write(
        "# using mapped file %s and %s as reference file\n" % (opt.filename, opt.ref)
    )
    if opt.minqual != 0:
        out.write(
            "# Quality filtering of bases with a Phred score < %d\n" % opt.minqual
        )
    out.write("# Std: strand of reads\n")
    out.write("Std\tKind\tLength\tOccurences \n")
    for (pe_or_se, strand), lengths in sorted(tab.items()):
        for length, count in sorted(lengths.items()):
            out.write("%s\t%s\t%d\t%d\n" % (strand, pe_or_se, length, count))


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
            return True

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
                return True

    if total < 0.01:
        logger.warning(
            "DNA damage levels are too low, the Bayesian computation should not be performed (%f < 0.01)",
            total,
        )

    return False


def _print_freq_table(table, columns, opt, out, offset=0):
    out.write("# table produced by mapDamage version %s\n" % __version__)
    out.write(
        "# using mapped file %s and %s as reference file\n" % (opt.filename, opt.ref)
    )
    if opt.minqual != 0:
        out.write("# Bases with a Phred score < %d were filtered out\n" % opt.minqual)
    out.write(
        "# Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads\n"
    )
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
