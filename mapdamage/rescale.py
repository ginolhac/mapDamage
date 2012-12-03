import csv
import sys
import os
import mapdamage
import pysam
import itertools
import math


def get_corr_prob(folder):
    """
    Reads the damage probability correction file, returns a 
    dictionary with this structure 
    position (one based)  -  CT  -  probability
                          -  GA  -  probability
    """
    try:
        fi_handle = csv.DictReader(open(os.path.join(folder,"Stats_out_MCMC_correct_prob.csv")))
        corr_prob = {}
        for line in fi_handle:
            if (corr_prob.has_key(line["Position"])):
                sys.exit('This file has multiple position definitions %s, line %d: %s' % \
                    (folder, fi_handle.line_num, corr_prob[line["Position"]]))
            else:
                corr_prob[int(line["Position"])] = {'C.T':float(line["C.T"]), 'G.A':float(line["G.A"])}
        return corr_prob
    except csv.Error as e:
        sys.exit('File %s, line %d: %s' % (os.path.join(folder,"Stats_out_MCMC_correct_prob.csv"), \
            fi_handle.line_num, e))


def corr_this_base(corr_prob, nt_seq, nt_ref, pos, length):
    """
    The position specific damaging correction, using the input 
    corr_prob dictionary holding the damage correcting values 
    nt_seq nucleotide in the sequence 
    nt_ref nucleotide in the reference
    pos relative position from the 5' end 
    length length of the sequence
    returns the correction probability for this particular set
    """
    if (pos == 0):
        # not using 0 based indexing
        raise SystemError
    if ( nt_seq == "T" and nt_ref == "C" ):
        # an C to T transition
        subs = "C.T"
    elif( nt_seq == "A" and nt_ref == "G" ):
        # an G to A transition
        subs = "G.A"
    else:
        # other transitions/transversions are not affected by damage
        return 0
    
    back_pos = pos-length-1 
    # position from 3' end

    if corr_prob.has_key(pos):
        p5_corr = corr_prob[pos][subs]
        # correction from 5' end
    else:
        p5_corr = 0

    if corr_prob.has_key(back_pos):
        p3_corr = corr_prob[back_pos][subs]
        # correction from 3' end
    else:
        p3_corr = 0

    if p5_corr != 0 and p3_corr != 0:
        # take a mean of the corrections this happens 
        # when the read is shorter than the length 
        # used for the statistical estimation.
        return (p5_corr+p3_corr)/2
    else:
        return max(p5_corr, p3_corr)


def rescale_qual_read(bam, read, ref, corr_prob, debug = False):
    """
    bam              a pysam bam object
    read             a pysam read object
    ref              a pysam fasta ref file
    reflengths       a dictionary holding the length of the references 
    corr_prob dictionary from get_corr_prob
    options            options from the command line parsing
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
    (seq, qual, refseq) = mapdamage.align.align_with_qual(read.cigar, \
        raw_seq, read.qqual, -100, refseq)
    length_read = len(raw_seq)
    length_align = len(seq)
    # reverse complement read and reference when mapped reverse strand
    if read.is_reverse:
        refseq = mapdamage.seq.revcomp(refseq)
        seq = mapdamage.seq.revcomp(seq)
        qual = qual[::-1]
    new_qual = [-100]*length_read
    pos_on_read = 0
    for (i, nt_seq, nt_ref, nt_qual) in itertools.izip(xrange(length_align), seq, refseq, qual):
        # rescale the quality according to the triplet position, 
        # pair of the reference and the sequence
        pdam = 1 - corr_this_base(corr_prob, nt_seq, nt_ref, pos_on_read + 1, length_read)
        pseq = 1 - 10**(-(float(ord(nt_qual))-float(33))/10)
        newp = pdam*pseq #This could be numerically unstable
        new_qual[pos_on_read] = chr(int(round(-10*math.log10(abs(1-newp)))+33))
        if nt_seq != "-":
            pos_on_read += 1
    # done with the aligned portion of the read 
    new_qual = "".join(new_qual)
    if read.is_reverse:
        new_qual = new_qual[::-1]

    if (read.cigar[0][0] == 4):
        # check for soft clipping at forward end 
        new_qual = read.qual[0:read.cigar[0][1]] + new_qual
    if (read.cigar[-1][0] == 4):
        # the same backwards
        new_qual = new_qual + read.qual[-read.cigar[-1][1]:]

    if debug:
        print ""
        print "ref-"+refseq 
        print "seq-"+seq
        print ""
        if (refseq != seq):
            print "Reference and the sequence are not the same"
            print [i for i, (left, right) in enumerate(zip(refseq, seq)) if left != right]
            print ""
        if read.is_reverse:
            print "rea-"+read.qual[::-1]
        else:
            print "rea-"+read.qual
        if read.is_reverse:
            print "new-"+new_qual[::-1]
        else:
            print "new-"+new_qual
        if (read.qual != new_qual):
            print "New and old qual are not the same"
            print [i for i, (left, right) in enumerate(zip(refseq, seq)) if left != right]
        print ""
        print ""
    # done
    read.qual = new_qual 

    return read


def rescale_qual(ref, options):
    """    
    ref                a pysam fasta ref file
    bam_filename       name of a BAM/SAM file to read
    fi                 file containing the csv with correction probabilities
    reflengths         dictionary with the reference lengths
    options            options from the command line parsing

    Iterates through BAM file, makes a new BAM file with rescaled qualities.
    """
    # open SAM/BAM file
    bam = pysam.Samfile(options.filename)

    # format the output filename
    file_name = os.path.splitext(options.filename)[0]
    out_file_name = os.path.basename(file_name)+"_rescaled.bam"

    bam_out = pysam.Samfile(os.path.join(options.folder, out_file_name), "wb", template = bam)
    
    corr_prob = get_corr_prob(options.folder)
    
    for hit in bam:
        hit = rescale_qual_read(bam, hit, ref, corr_prob)
        if hit.is_paired:
            sys.stderr.write("Cannot rescale paired end reads in this versio\n")
            sys.exit(1)
        bam_out.write(hit)
    bam.close()
    bam_out.close()
    if not options.quiet:
        print("Done with rescaling.")
    if options.verbose:
        print("Rescaled BAM: %s" % os.path.join(options.folder, out_file_name))

