#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Copyright (c) 2012  Aurélien Ginolhac, Mikkel Schubert, Hákon Jónsson
and Ludovic Orlando

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.

plot and quantify damage patterns from a SAM/BAM file

:Authors: Aurélien Ginolhac, Mikkel Schubert, Hákon Jónsson, Ludovic Orlando
:Date: November 2012
:Type: tool
:Input: SAM/BAM
:Output: tabulated tables, pdf
"""
import logging
import time
import sys
import os

import coloredlogs
import pysam

import mapdamage

_LOG_FORMAT = "%(asctime)s %(name)s %(levelname)s %(message)s"


def main(argv):
    start_time = time.time()

    coloredlogs.install(fmt=_LOG_FORMAT)
    logger = logging.getLogger(__name__)

    options = mapdamage.parseoptions.options(argv)
    if options is None:
        logging.error("Option parsing failed, terminating the program")
        return 1

    handler = logging.FileHandler(os.path.join(options.folder, "Runtime_log.txt"))
    formatter = logging.Formatter(_LOG_FORMAT)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    logging.getLogger().setLevel(options.log_level)
    for handler in logging.getLogger().handlers:
        handler.setLevel(options.log_level)

    logger.info("Started with the command: " + " ".join(sys.argv))

    # plot using R if results folder already done
    if options.plot_only:
        if options.no_r:
            logger.error("Cannot use plot damage patterns if R is missing, terminating")
            return 1
        else:
            if not mapdamage.rscript.misincorporation_plot(options):
                return 1

            if not mapdamage.rscript.length_distribution_plot(options):
                return 1

            return 0

    # run the Bayesian estimation if the matrix construction is done
    if options.stats_only:
        # does not work for very low damage levels
        if mapdamage.tables.check_table_and_warn_if_dmg_freq_is_low(options.folder):
            logger.error("Cannot use the Bayesian estimation, terminating the program")
            return 1
        else:
            # before running the Bayesian estimation get the base composition
            path_to_basecomp = os.path.join(options.folder, "dnacomp_genome.csv")
            if os.path.isfile(path_to_basecomp):
                # Try to read the base composition file
                mapdamage.composition.read_base_comp(path_to_basecomp)
            else:
                # Construct the base composition file
                mapdamage.composition.write_base_comp(options.ref, path_to_basecomp)

            if not mapdamage.rscript.perform_bayesian_estimates(options):
                return 1

            return 0

    # fetch all references and associated lengths in nucleotides
    try:
        ref = pysam.Fastafile(options.ref)
    except IOError as error:
        logger.error("Could not open the reference file %r: %e", options.ref, error)
        raise

    # rescale the qualities
    if options.rescale_only:
        logger.info("Starting rescaling...")
        mapdamage.rescale.rescale_qual(ref, options)
        return 0

    # open SAM/BAM file
    reader = mapdamage.reader.BAMReader(
        filepath=options.filename,
        downsample_to=options.downsample,
        downsample_seed=options.downsample_seed,
    )

    if reader.is_stream and options.rescale:
        # rescaling is not possible on a streasm, since we need to read it twice
        logger.warning("Reading from stream, rescaling is disabled")
        options.rescale = False

    reflengths = reader.get_references()
    # check if references in SAM/BAM are the same in the fasta reference file
    fai_lengths = mapdamage.seq.read_fasta_index(options.ref + ".fai")
    if not fai_lengths:
        return 1
    elif not mapdamage.seq.compare_sequence_dicts(fai_lengths, reflengths):
        return 1
    elif (len(reflengths) >= 1000) and not options.merge_reference_sequences:
        logger.warning(
            "Alignment contains a large number of reference sequences (%i)!",
            len(reflengths),
        )
        logger.warning("This may lead to excessive memory/disk usage.")
        logger.warning("Consider using --merge-reference-sequences")

    refnames = reader.handle.references
    if options.merge_reference_sequences:
        refnames = ["*"]

    # for misincorporation patterns, record mismatches
    misincorp = mapdamage.tables.initialize_mut(refnames, options.length)
    # for fragmentation patterns, record base compositions
    dnacomp = mapdamage.tables.initialize_comp(refnames, options.around, options.length)
    # for length distributions
    lgdistrib = mapdamage.tables.initialize_lg()

    logger.info("Reading from '%s'", options.filename)
    if options.minqual != 0:
        logger.info("Filtering out bases with a Phred score < %d", options.minqual)
    logger.debug(
        "%d references are assumed in SAM/BAM file, for a total of %d nucleotides",
        len(reflengths),
        sum(reflengths.values()),
    )
    logger.info("Writing results to '%s/'", options.folder)

    # main loop
    counter = 0
    warned_about_quals = False
    for read in reader:
        counter += 1

        # external coordinates 5' and 3' , 3' is 1-based offset
        coordinate = mapdamage.align.get_coordinates(read)
        # record aligned length for single-end reads
        mapdamage.seq.record_length(read, coordinate, lgdistrib)
        # fetch reference name, chromosome or contig names
        chrom = reader.handle.getrname(read.tid)

        (before, after) = mapdamage.align.get_around(
            coordinate, chrom, reflengths, options.around, ref
        )
        refseq = ref.fetch(chrom, min(coordinate), max(coordinate)).upper()
        # read.query contains aligned sequences while read.seq is the read itself
        seq = read.query

        # add gaps according to the cigar string, do it for qualities if filtering options is on
        if not (options.minqual and read.qual):
            if options.minqual and not warned_about_quals:
                logger.warning(
                    "Reads without PHRED scores found; cannot filter by --min-basequal"
                )
                warned_about_quals = True

            (seq, refseq) = mapdamage.align.align(read.cigar, seq, refseq)
        else:
            # add gaps to qualities and mask read and reference nucleotides if below desired threshold
            (seq, _, refseq) = mapdamage.align.align_with_qual(
                read.cigar, seq, read.qqual, options.minqual, refseq
            )

        # reverse complement read and reference when mapped reverse strand
        if read.is_reverse:
            refseq = mapdamage.seq.revcomp(refseq)
            seq = mapdamage.seq.revcomp(seq)
            beforerev = mapdamage.seq.revcomp(after)
            after = mapdamage.seq.revcomp(before)
            before = beforerev

        if options.merge_reference_sequences:
            chrom = "*"

        # record soft clipping when present
        mapdamage.align.record_soft_clipping(read, misincorp[chrom], options.length)

        # count misincorparations by comparing read and reference base by base
        mapdamage.align.get_mis(
            read, seq, refseq, chrom, options.length, misincorp, "5p"
        )
        # do the same with sequences align to 3'-ends
        mapdamage.align.get_mis(
            read, seq[::-1], refseq[::-1], chrom, options.length, misincorp, "3p"
        )
        # compute base composition for reads
        mapdamage.composition.count_read_comp(read, chrom, options.length, dnacomp)

        # compute base composition for genomic regions
        mapdamage.composition.count_ref_comp(read, chrom, before, after, dnacomp)

        if counter % 50000 == 0:
            logger.debug("%10d filtered alignments processed", counter)

    logger.debug("Done. %d filtered alignments processed", counter)
    logger.debug("BAM read in %f seconds", time.time() - start_time)

    # close file handles
    reader.close()

    # output results, write summary tables to disk
    with open(os.path.join(options.folder, "misincorporation.txt"), "w") as fmut:
        mapdamage.tables.print_mut(misincorp, options, fmut)
    with open(os.path.join(options.folder, "dnacomp.txt"), "w") as fcomp:
        mapdamage.tables.print_comp(dnacomp, options, fcomp)
    with open(os.path.join(options.folder, "lgdistribution.txt"), "w") as flg:
        mapdamage.tables.print_lg(lgdistrib, options, flg)

    # plot using R
    if not options.no_r:
        if not mapdamage.rscript.misincorporation_plot(options):
            return 1

        if not mapdamage.rscript.length_distribution_plot(options):
            return 1

    # raises a warning for very low damage levels
    if mapdamage.tables.check_table_and_warn_if_dmg_freq_is_low(options.folder):
        options.no_stats = True

    # run the Bayesian estimation
    if not options.no_stats:
        # before running the Bayesian estimation get the base composition
        mapdamage.composition.write_base_comp(
            options.ref, os.path.join(options.folder, "dnacomp_genome.csv")
        )

        if not mapdamage.rscript.perform_bayesian_estimates(options):
            return 1

    # rescale the qualities
    if options.rescale:
        mapdamage.rescale.rescale_qual(ref, options)

    # need the fasta reference still open for rescaling
    ref.close()

    # log the time it took
    logger.info("Successful run")
    logger.debug("Run completed in %f seconds" % (time.time() - start_time,))

    return 0


def entry_point():
    return main(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
