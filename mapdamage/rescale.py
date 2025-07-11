import csv
import logging
import math
import os
import sys
import time

import pysam

import mapdamage


def phred_pval_to_char(pval):
    """Transforming error rate to ASCII character using the Phred scale"""
    return chr(int(round(-10 * math.log10(abs(pval))) + 33))


def phred_char_to_pval(ch):
    """Transforming ASCII character in the Phred scale to the error rate"""
    return 10 ** (-(float(ord(ch)) - float(33)) / 10)


def get_corr_prob(folder, rescale_length_5p, rescale_length_3p):
    """
    Reads the damage probability correction file, returns a
    dictionary with this structure
    position (one based)  -  CT  -  probability
                          -  GA  -  probability
    """
    full_path = os.path.join(folder, "Stats_out_MCMC_correct_prob.csv")
    if not os.path.isfile(full_path):
        sys.exit(
            "Missing file, the file \n\tStats_out_MCMC_correct_prob.csv\nshould be in the folder\n\t"
            + folder
            + "\nDid you run the MCMC estimates of the parameters?"
        )
    try:
        with open(full_path) as fi:
            fi_handle = csv.DictReader(fi)
            corr_prob = {}
            for line in fi_handle:
                if line["Position"] in corr_prob:
                    sys.exit(
                        "This file has multiple position definitions %s, line %d: %s"
                        % (folder, fi_handle.line_num, corr_prob[line["Position"]])
                    )
                else:
                    corr_prob[int(line["Position"])] = {
                        "C.T": float(line["C.T"]),
                        "G.A": float(line["G.A"]),
                    }

            # Exclude probabilities for positions outside of user-specified region
            for key in list(corr_prob.keys()):
                if key < -rescale_length_3p or key > rescale_length_5p:
                    corr_prob.pop(key)
            return corr_prob
    except csv.Error as e:
        sys.exit(
            "File %s, line %d: %s"
            % (
                os.path.join(folder, "Stats_out_MCMC_correct_prob.csv"),
                fi_handle.line_num,
                e,
            )
        )


def corr_this_base(corr_prob, nt_seq, nt_ref, pos, length, direction="both"):
    """
    The position specific damaging correction, using the input
    corr_prob dictionary holding the damage correcting values
    nt_seq nucleotide in the sequence
    nt_ref nucleotide in the reference
    pos relative position from the 5' end
    length length of the sequence
    direction which end to consider the rescaling
    returns the correction probability for this particular set
    """
    if pos == 0:
        # not using 0 based indexing
        raise SystemError
    if nt_seq == "T" and nt_ref == "C":
        # an C to T transition
        subs = "C.T"
    elif nt_seq == "A" and nt_ref == "G":
        # an G to A transition
        subs = "G.A"
    else:
        # other transitions/transversions are not affected by damage
        return 0

    back_pos = pos - length - 1
    # position from 3' end

    if pos in corr_prob:
        p5_corr = corr_prob[pos][subs]
        # correction from 5' end
    else:
        p5_corr = 0

    if back_pos in corr_prob:
        p3_corr = corr_prob[back_pos][subs]
        # correction from 3' end
    else:
        p3_corr = 0

    if direction == "forward":
        return p5_corr
    elif direction == "backward":
        return p3_corr
    elif direction == "both":
        if pos < abs(back_pos):
            # then we use the forward correction
            return p5_corr
        else:
            # else the backward correction
            return p3_corr
    else:
        # this should not happen
        raise SystemExit("Abnormal direction in the rescaling procedure")


def initialize_subs():
    """Initialize a substitution table, to track the expected substitution counts"""
    per_qual = dict(list(zip(list(range(0, 130)), [0] * 130)))
    subs = {
        "CT-before": per_qual.copy(),
        "TC-before": per_qual.copy(),
        "GA-before": per_qual.copy(),
        "AG-before": per_qual.copy(),
        "CT-after": per_qual.copy(),
        "TC-after": per_qual.copy(),
        "GA-after": per_qual.copy(),
        "AG-after": per_qual.copy(),
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
        "CT-pvals": 0.0,
        "CT-pvals_before": 0.0,
        "TC-pvals": 0.0,
        "GA-pvals": 0.0,
        "GA-pvals_before": 0.0,
        "AG-pvals": 0.0,
    }
    return subs


def record_subs(subs, nt_seq, nt_ref, nt_qual, nt_newqual, prob_corr):
    """record the expected substitution change, prob_corr is the exact version for nt_qual"""
    if nt_seq == "T" and nt_ref == "C":
        sub_type = "CT"
        subs["CT-pvals"] += prob_corr
        subs["CT-pvals_before"] += 1 - phred_char_to_pval(nt_qual)
    elif nt_seq == "A" and nt_ref == "G":
        sub_type = "GA"
        subs["GA-pvals"] += prob_corr
        subs["GA-pvals_before"] += 1 - phred_char_to_pval(nt_qual)
    elif nt_seq == "C" and nt_ref == "T":
        sub_type = "TC"
        subs["TC-pvals"] += 1 - phred_char_to_pval(nt_qual)
        if nt_qual != nt_newqual:
            raise SystemError(
                "Internal error: rescaling qualities for the wrong transitions"
            )
    elif nt_seq == "G" and nt_ref == "A":
        sub_type = "AG"
        subs["AG-pvals"] += 1 - phred_char_to_pval(nt_qual)
        if nt_qual != nt_newqual:
            raise SystemError(
                "Internal error: rescaling qualities for the wrong transitions"
            )
    else:
        sub_type = "NN"
    if sub_type != "NN":
        # record only transitions
        subs[sub_type + "-before"][int(ord(nt_qual)) - 33] += 1
        subs[sub_type + "-after"][int(ord(nt_newqual)) - 33] += 1
    if nt_ref in ["A", "C", "G", "T"]:
        subs[nt_ref] += 1


def qual_summary_subs(subs):
    """Calculates summary statistics for the substitution table subs"""
    for i in [
        "CT-before",
        "TC-before",
        "GA-before",
        "AG-before",
        "CT-after",
        "TC-after",
        "GA-after",
        "AG-after",
    ]:
        for lv in [0, 10, 20, 30, 40]:
            for qv in subs[i]:
                if qv >= lv:
                    key = i + "-Q" + str(lv)
                    if key in subs:
                        subs[key] += subs[i][qv]
                    else:
                        subs[key] = subs[i][qv]


def print_subs(subs):
    """Print the substitution table"""
    print(
        "\tThe expected substitution frequencies before and after scaling using the scaled qualities as probalities:"
    )
    if subs["C"] != 0:
        # the special case of no substitutions
        print(
            "\tCT\t"
            + str(subs["CT-pvals_before"] / subs["C"])
            + "\t\t"
            + str(subs["CT-pvals"] / subs["C"])
        )
    else:
        print("\tCT\tNA\t\tNA")
    if subs["T"] != 0:
        print(
            "\tTC\t"
            + str(subs["TC-pvals"] / subs["T"])
            + "\t\t"
            + str(subs["TC-pvals"] / subs["T"])
        )
    else:
        print("\tTC\tNA\t\tNA")
    if subs["G"] != 0:
        print(
            "\tGA\t"
            + str(subs["GA-pvals_before"] / subs["G"])
            + "\t\t"
            + str(subs["GA-pvals"] / subs["G"])
        )
    else:
        print("\tGA\tNA\t\tNA")
    if subs["A"] != 0:
        print(
            "\tAG\t"
            + str(subs["AG-pvals"] / subs["A"])
            + "\t\t"
            + str(subs["AG-pvals"] / subs["A"])
        )
    else:
        print("\tAG\tNA\t\tNA")
    print("\tQuality metrics before and after scaling")
    print("\tCT-Q0 \t" + str(subs["CT-before-Q0"]) + "\t\t" + str(subs["CT-after-Q0"]))
    print(
        "\tCT-Q10 \t" + str(subs["CT-before-Q10"]) + "\t\t" + str(subs["CT-after-Q10"])
    )
    print(
        "\tCT-Q20 \t" + str(subs["CT-before-Q20"]) + "\t\t" + str(subs["CT-after-Q20"])
    )
    print(
        "\tCT-Q30 \t" + str(subs["CT-before-Q30"]) + "\t\t" + str(subs["CT-after-Q30"])
    )
    print(
        "\tCT-Q40 \t" + str(subs["CT-before-Q40"]) + "\t\t" + str(subs["CT-after-Q40"])
    )
    print("\tGA-Q0 \t" + str(subs["GA-before-Q0"]) + "\t\t" + str(subs["GA-after-Q0"]))
    print(
        "\tGA-Q10 \t" + str(subs["GA-before-Q10"]) + "\t\t" + str(subs["GA-after-Q10"])
    )
    print(
        "\tGA-Q20 \t" + str(subs["GA-before-Q20"]) + "\t\t" + str(subs["GA-after-Q20"])
    )
    print(
        "\tGA-Q30 \t" + str(subs["GA-before-Q30"]) + "\t\t" + str(subs["GA-after-Q30"])
    )
    print(
        "\tGA-Q40 \t" + str(subs["GA-before-Q40"]) + "\t\t" + str(subs["GA-after-Q40"])
    )


def rescale_qual_read(bam, read, ref, corr_prob, subs, debug=False, direction="both"):
    """
    bam              a pysam bam object
    read             a pysam read object
    ref              a pysam fasta ref file
    reflengths       a dictionary holding the length of the references
    subs             a dictionary holding the corrected number of substitution before and after scaling
    corr_prob dictionary from get_corr_prob
    returns a read with rescaled quality score

    Iterates through the read and reference, rescales the quality
    according to corr_prob
    """
    if not debug:
        # no need to log when unit testing
        logger = logging.getLogger(__name__)
    raw_seq = read.query
    # external coordinates 5' and 3' , 0-based offset
    coordinate = mapdamage.align.get_coordinates(read)
    # fetch reference name, chromosome or contig names
    chrom = bam.getrname(read.tid)
    refseq = ref.fetch(chrom, min(coordinate), max(coordinate)).upper()
    # add gaps to qualities and mask read and reference nucleotides if below desired threshold
    (seq, qual, refseq) = mapdamage.align.align_with_qual(
        read.cigar, raw_seq, read.qqual, -100, refseq
    )
    length_read = len(raw_seq)
    length_align = len(seq)
    # reverse complement read and reference when mapped reverse strand
    if read.is_reverse:
        refseq = mapdamage.seq.revcomp(refseq)
        seq = mapdamage.seq.revcomp(seq)
        qual = qual[::-1]
    new_qual = [-100] * length_read
    pos_on_read = 0
    number_of_rescaled_bases = 0.0
    for i, nt_seq, nt_ref, nt_qual in zip(range(length_align), seq, refseq, qual):
        # rescale the quality according to the triplet position,
        # pair of the reference and the sequence
        if (nt_seq == "T" and nt_ref == "C") or (nt_seq == "A" and nt_ref == "G"):
            # need to rescale this subs.
            pdam = 1 - corr_this_base(
                corr_prob,
                nt_seq,
                nt_ref,
                pos_on_read + 1,
                length_read,
                direction=direction,
            )
            pseq = 1 - phred_char_to_pval(nt_qual)
            newp = pdam * pseq  # this could be numerically unstable
            newq = phred_pval_to_char(1 - newp)
            number_of_rescaled_bases += 1 - pdam
        else:
            # don't rescale, other bases
            newp = 1 - phred_char_to_pval(nt_qual)
            newq = nt_qual
        if pos_on_read < length_read:
            new_qual[pos_on_read] = newq
            record_subs(subs, nt_seq, nt_ref, nt_qual, new_qual[pos_on_read], newp)
            if nt_seq != "-":
                pos_on_read += 1
            # done with the aligned portion of the read
        else:
            if not debug:
                logger.warning(
                    "Warning: The alignment of the read is longer than the actual read %s",
                    (read.qname),
                )
            break
    new_qual = "".join(new_qual)

    if read.is_reverse:
        new_qual = new_qual[::-1]
    if read.cigar[0][0] == 4:
        # check for soft clipping at forward end
        new_qual = read.qual[0 : read.cigar[0][1]] + new_qual
    if read.cigar[-1][0] == 4:
        # the same backwards
        new_qual = new_qual + read.qual[-read.cigar[-1][1] :]

    read.qual = new_qual
    # truncate this to 5 digits
    number_of_rescaled_bases = float("%.5f" % number_of_rescaled_bases)

    if read.has_tag("MR"):
        raise SystemExit("Read: %s already has a MR tag, can't rescale" % read)

    read.set_tag("MR", number_of_rescaled_bases, "f")

    return read


def rescale_qual(ref, options, debug=False):
    """
    ref                a pysam fasta ref file
    bam_filename       name of a BAM/SAM file to read
    fi                 file containing the csv with correction probabilities
    reflengths         dictionary with the reference lengths
    options            options from the command line parsing

    Iterates through BAM file, makes a new BAM file with rescaled qualities.
    """
    if not debug:
        # no need to log when unit testing
        logger = logging.getLogger(__name__)
        logger.info(
            "Rescaling BAM: '%s' -> '%s'" % (options.filename, options.rescale_out)
        )
    start_time = time.time()

    # open SAM/BAM files
    bam = pysam.Samfile(options.filename)
    if debug:
        write_mode = "wh"
    else:
        write_mode = "wb"
    bam_out = pysam.Samfile(options.rescale_out, write_mode, template=bam)
    corr_prob = get_corr_prob(
        options.folder,
        rescale_length_5p=options.rescale_length_5p,
        rescale_length_3p=options.rescale_length_3p,
    )
    subs = initialize_subs()
    first_pair = True
    number_of_non_proper_pairs = 0

    for hit in bam:
        if hit.is_unmapped:
            pass
        elif not hit.qual and not debug:
            logger.warning(
                "Cannot rescale base PHRED scores for read '%s'; no scores assigned."
                % hit.qname
            )
        elif hit.is_paired:
            if first_pair and not debug:
                # assuming the ends are non-overlapping
                logger.warning(
                    "Warning! Assuming the pairs are non-overlapping, facing inwards and correctly paired."
                )
                first_pair = False
            # 5p --------------> 3p
            # 3p <-------------- 5p
            # pair 1 (inwards)
            # 5p ---->
            #             <---- 5p
            #     A         B
            # pair 2 (outwards), this happens if the reference is RC this is not supported
            #             ----> 3p
            # 3p <----
            #     A         B
            # Correct outwards pairs from the 3p and inwards pairs with the 5p end
            if (
                (not hit.is_reverse)
                and hit.mate_is_reverse
                and (hit.pnext > hit.pos)
                and hit.tid == hit.mrnm
            ):
                # the inwards case mate A
                hit = rescale_qual_read(
                    bam, hit, ref, corr_prob, subs, direction="forward", debug=debug
                )
            elif (
                hit.is_reverse
                and (not hit.mate_is_reverse)
                and (hit.pnext < hit.pos)
                and hit.tid == hit.mrnm
            ):
                # the inwards case mate B
                hit = rescale_qual_read(
                    bam, hit, ref, corr_prob, subs, direction="forward", debug=debug
                )
            else:
                number_of_non_proper_pairs += 1
                # cannot do much with conflicting pairing information
        else:
            hit = rescale_qual_read(bam, hit, ref, corr_prob, subs, debug=debug)

        bam_out.write(hit)
    if number_of_non_proper_pairs != 0 and not debug:
        logger.warning(
            "Number of non-rescaled reads due to improper pairing:  %d"
            % number_of_non_proper_pairs
        )
    if subs["TC-before"] != subs["TC-after"] or subs["AG-before"] != subs["AG-after"]:
        sys.exit(
            "Qualities for T.C and A.G transitions should not change in the rescaling, please contact the authors."
        )
    qual_summary_subs(subs)
    bam.close()
    bam_out.close()
    if not options.quiet:
        print_subs(subs)
    if not debug:
        logger.debug("Rescaling completed in %f seconds" % (time.time() - start_time,))
