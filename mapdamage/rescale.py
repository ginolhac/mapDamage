import csv
import mapdamage
import pysam
import math
import logging
import time


class RescaleError(RuntimeError):
    pass


def _phred_pval_to_char(pval):
    """ Transforming error rate to ASCII character using the Phred scale"""
    return chr(int(round(-10 * math.log10(abs(pval))) + 33))


def _phred_char_to_pval(ch):
    """ Transforming ASCII character in the Phred scale to the error rate"""
    return 10 ** (-(float(ord(ch)) - float(33)) / 10)


def _get_corr_prob(filepath, rescale_length_5p, rescale_length_3p):
    """ Reads the damage probability correction table, and returns a dictionary with the
    structure {(ref_nt, read_nt, position): probability}
    """
    logger = logging.getLogger(__name__)
    logger.info("Reading corrected probabilities from '%s'", filepath)

    try:
        with filepath.open(newline="") as handle:
            reader = csv.DictReader(handle, strict=True)
            corr_prob = {}
            for line in reader:
                position = int(line["Position"])

                # Exclude probabilities for positions outside of user-specified region
                if -rescale_length_3p <= position <= rescale_length_5p:
                    corr_prob[("C", "T", position)] = float(line["C.T"])
                    corr_prob[("G", "A", position)] = float(line["G.A"])

            return corr_prob
    except FileNotFoundError:
        raise RescaleError("File does not exist; please re-run mapDamage")
    except csv.Error as error:
        raise RescaleError("Error while reading line %d: %s" % (reader.line_num, error))


def _corr_this_base(corr_prob, nt_seq, nt_ref, pos, length, direction="both"):
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

    # position from 3' end
    back_pos = pos - length - 1

    if direction == "both":
        if pos >= abs(back_pos):
            pos = back_pos
    elif direction == "reverse":
        pos = back_pos
    elif direction != "forward":
        # this should not happen
        raise RescaleError(
            "Abnormal direction in the rescaling procedure (%r); please submit a bug-"
            "report on github" % (direction,)
        )

    return corr_prob.get((nt_ref, nt_seq, pos), 0)


def _initialize_subs():
    """Initialize a substitution table, to track the expected substitution counts"""
    per_qual = dict(zip(range(130), [0] * 130))
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


def _record_subs(subs, nt_seq, nt_ref, nt_qual, nt_newqual, prob_corr):
    """ record the expected substitution change, prob_corr is the excact version for nt_qual"""
    if nt_seq == "T" and nt_ref == "C":
        sub_type = "CT"
        subs["CT-pvals"] += prob_corr
        subs["CT-pvals_before"] += 1 - _phred_char_to_pval(nt_qual)
    elif nt_seq == "A" and nt_ref == "G":
        sub_type = "GA"
        subs["GA-pvals"] += prob_corr
        subs["GA-pvals_before"] += 1 - _phred_char_to_pval(nt_qual)
    elif nt_seq == "C" and nt_ref == "T":
        sub_type = "TC"
        subs["TC-pvals"] += 1 - _phred_char_to_pval(nt_qual)
        if nt_qual != nt_newqual:
            raise SystemError(
                "Internal error: rescaling qualities for the wrong transitions"
            )
    elif nt_seq == "G" and nt_ref == "A":
        sub_type = "AG"
        subs["AG-pvals"] += 1 - _phred_char_to_pval(nt_qual)
        if nt_qual != nt_newqual:
            raise SystemError(
                "Internal error: rescaling qualities for the wrong transitions"
            )
    else:
        sub_type = "NN"
    if sub_type != "NN":
        # record only transitions
        subs[sub_type + "-before"][ord(nt_qual) - 33] += 1
        subs[sub_type + "-after"][ord(nt_newqual) - 33] += 1
    if nt_ref in ["A", "C", "G", "T"]:
        subs[nt_ref] += 1


def _qual_summary_subs(subs):
    """Calculates summary statistics for the substition table subs"""
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


def _print_subs(subs):
    """Print the substition table"""
    log = logging.getLogger(__name__).info
    log("Expected substition frequencies before and after rescaling:")
    for sub in ("CT", "TC", "GA", "AG"):
        base_count = subs[sub[0]]

        if base_count:
            pvals_key = sub + "-pvals"
            pvals = subs[pvals_key]
            pvals_before = subs.get(pvals_key + "_before", pvals)

            log(
                "    %s>%s    %.4f    %.4f",
                sub[0],
                sub[1],
                pvals_before / base_count,
                pvals / base_count,
            )
        else:
            log("\t%s\tNA\t\tNA", sub)

    log("Quality metrics before and after scaling:")
    for sub in ("CT", "GA"):
        for qual in (0, 10, 20, 30, 40):
            before = subs["%s-before-Q%i" % (sub, qual)]
            after = subs["%s-after-Q%i" % (sub, qual)]

            log("    %s-Q%02i% 10i% 10i", sub, qual, before, after)


def _rescale_qual_read(bam, read, ref, corr_prob, subs, direction="both"):
    """
    bam              a pysam bam object
    read             a pysam read object
    ref              a pysam fasta ref file
    reflengths       a dictionary holding the length of the references
    subs             a dictionary holding the corrected number of substition before and after scaling
    corr_prob dictionary from _get_corr_prob
    returns a read with rescaled quality score

    Iterates through the read and reference, rescales the quality
    according to corr_prob
    """
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
    for (_, nt_seq, nt_ref, nt_qual) in zip(range(length_align), seq, refseq, qual):
        # rescale the quality according to the triplet position,
        # pair of the reference and the sequence
        if (nt_seq == "T" and nt_ref == "C") or (nt_seq == "A" and nt_ref == "G"):
            # need to rescale this subs.
            pdam = 1 - _corr_this_base(
                corr_prob,
                nt_seq,
                nt_ref,
                pos_on_read + 1,
                length_read,
                direction=direction,
            )
            pseq = 1 - _phred_char_to_pval(nt_qual)
            newp = pdam * pseq  # this could be numerically unstable
            newq = _phred_pval_to_char(1 - newp)
            number_of_rescaled_bases += 1 - pdam
        else:
            # don't rescale, other bases
            newp = 1 - _phred_char_to_pval(nt_qual)
            newq = nt_qual
        if pos_on_read < length_read:
            new_qual[pos_on_read] = newq
            _record_subs(subs, nt_seq, nt_ref, nt_qual, new_qual[pos_on_read], newp)
            if nt_seq != "-":
                pos_on_read += 1
            # done with the aligned portion of the read
        else:
            logger = logging.getLogger(__name__)
            logger.warning(
                "The aligment of the read is longer than the actual read %s",
                read.qname,
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


def _rescale_qual_core(ref, options):
    """Iterates through BAM file, writing new BAM file with rescaled qualities.
    """
    corr_prob = _get_corr_prob(
        filepath=options.folder / "Stats_out_MCMC_correct_prob.csv",
        rescale_length_5p=options.rescale_length_5p,
        rescale_length_3p=options.rescale_length_3p,
    )

    n_pairs = 0
    n_improper_pairs = 0
    n_reads_without_quals = 0
    subs = _initialize_subs()

    with pysam.AlignmentFile(options.filename) as bam_in:
        with pysam.AlignmentFile(options.rescale_out, "wb", template=bam_in) as bam_out:
            for hit in bam_in:
                if hit.is_unmapped:
                    pass
                elif not hit.qual:
                    n_reads_without_quals += 1
                elif hit.is_paired:
                    n_pairs += 1
                    # 5p --------------> 3p
                    # 3p <-------------- 5p
                    # pair 1 (inwards)
                    # 5p ---->
                    #             <---- 5p
                    #     A         B
                    # pair 2 (outwards); this is not supported
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
                        hit = _rescale_qual_read(
                            bam_in, hit, ref, corr_prob, subs, direction="forward"
                        )
                    elif (
                        hit.is_reverse
                        and (not hit.mate_is_reverse)
                        and (hit.pnext < hit.pos)
                        and hit.tid == hit.mrnm
                    ):
                        # the inwards case mate B
                        hit = _rescale_qual_read(
                            bam_in, hit, ref, corr_prob, subs, direction="forward"
                        )
                    else:
                        n_improper_pairs += 1
                        # cannot do much with conflicting pairing information
                else:
                    hit = _rescale_qual_read(bam_in, hit, ref, corr_prob, subs)

                bam_out.write(hit)

    logger = logging.getLogger(__name__)
    if n_pairs:
        logger.warning(
            "Processed %i paired reads, assumed to be non-overlapping, facing inwards "
            "and correctly paired; %i of these were excluded as improperly paired.",
            n_pairs,
            n_improper_pairs,
        )

    if n_reads_without_quals:
        logger.warning("Skipped %i reads without quality scores", n_reads_without_quals)

    if subs["TC-before"] != subs["TC-after"] or subs["AG-before"] != subs["AG-after"]:
        raise RescaleError(
            "Qualities for T.C and A.G transitions should not change in the rescaling. "
            "Please file a bug on github."
        )

    _qual_summary_subs(subs)
    _print_subs(subs)


def rescale_qual(ref, options):
    logger = logging.getLogger(__name__)
    logger.info("Rescaling BAM: '%s' -> '%s'", options.filename, options.rescale_out)
    start_time = time.time()

    try:
        _rescale_qual_core(ref, options)
    except RescaleError as error:
        logger.error("%s", error)
        return 1
    except Exception as error:
        logger.error("Unhandled exception: %s", error)
        return 1

    logger.debug("Rescaling completed in %f seconds", time.time() - start_time)
    return 0
