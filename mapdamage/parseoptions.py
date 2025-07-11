import os
import sys
from optparse import SUPPRESS_HELP, OptionGroup, OptionParser

from mapdamage.rscript import check_R_lib
from mapdamage.version import __version__


def file_exist(filename):
    if os.path.exists(filename) and not os.path.isdir(filename):
        return True
    elif filename == "-":
        return True
    else:
        sys.stderr.write("Error: '%s' is not a valid file\n" % (filename))
        return None


def whereis(program):
    for path in os.environ.get("PATH", "").split(":"):
        if os.path.exists(os.path.join(path, program)) and not os.path.isdir(
            os.path.join(path, program)
        ):
            return os.path.join(path, program)
    return None


def check_py_version():
    req_version = (2, 6)
    cur_version = sys.version_info

    if cur_version >= req_version:
        return True
    else:
        sys.stderr.write(
            "Your Python interpreter is too old."
            "Please consider upgrading to at least %d.%d\n"
            % (req_version[0], req_version[1])
        )
        return None


def options():
    parser = OptionParser(
        "%prog [options] -i BAMfile -r reference.fasta\n\nUse option -h or --help for help",
        version=__version__,
        epilog="report bugs to aginolhac@snm.ku.dk, MSchubert@snm.ku.dk or jonsson.hakon@gmail.com",
    )

    args = OptionGroup(parser, "Input files")
    args.add_option(
        "-i",
        "--input",
        help="SAM/BAM file, must contain a valid header, use '-' for reading a BAM from stdin",
        action="store",
        type="string",
        dest="filename",
    )
    args.add_option(
        "-r",
        "--reference",
        help="Reference file in FASTA format",
        action="store",
        dest="ref",
    )

    parser.add_option_group(args)
    group = OptionGroup(parser, "General options")
    group.add_option(
        "-n",
        "--downsample",
        help="Downsample to a randomly selected fraction of the reads (if 0 < DOWNSAMPLE < 1), or "
        "a fixed number of randomly selected reads (if DOWNSAMPLE >= 1). By default, no downsampling is performed.",
        type=float,
        default=None,
    )
    group.add_option(
        "--downsample-seed",
        help="Seed value to use for downsampling. See documentation for py module 'random' for default behavior.",
        type=int,
        default=None,
    )
    group.add_option(
        "--merge-reference-sequences",
        help="Ignore reference sequence names when tabulating reads (using '*' instead). "
        "Useful for alignments with a large number of reference sequences, which may otherwise result in excessive "
        "memory or disk usage due to the number of tables generated.",
        default=False,
        action="store_true",
    )
    group.add_option(
        "-l",
        "--length",
        dest="length",
        help="read length, in nucleotides to consider [%default]",
        type=int,
        default=70,
        action="store",
    )
    group.add_option(
        "-a",
        "--around",
        dest="around",
        help="nucleotides to retrieve before/after reads [%default]",
        type=int,
        default=10,
        action="store",
    )
    group.add_option(
        "-Q",
        "--min-basequal",
        dest="minqual",
        help="minimum base quality Phred score considered, Phred-33 assumed [%default]",
        type=int,
        default=0,
        action="store",
    )
    group.add_option(
        "-d",
        "--folder",
        help="folder name to store results [results_FILENAME]",
        action="store",
        type="string",
        dest="folder",
    )
    group.add_option(
        "-f",
        "--fasta",
        dest="fasta",
        help="Write alignments in a FASTA file",
        default=False,
        action="store_true",
    )
    group.add_option(
        "--plot-only",
        dest="plot_only",
        help="Run only plotting from a valid result folder",
        default=False,
        action="store_true",
    )
    group.add_option(
        "-q",
        "--quiet",
        dest="quiet",
        help="Disable any output to stdout",
        default=False,
        action="store_true",
    )
    group.add_option(
        "-v",
        "--verbose",
        dest="verbose",
        help="Display progression information during parsing",
        default=False,
        action="store_true",
    )
    group.add_option(
        "--no-plot", dest="no_r", help=SUPPRESS_HELP, default=False, action="store_true"
    )
    parser.add_option_group(group)

    # options for plotting damage patterns
    group2 = OptionGroup(parser, "Options for graphics")
    group2.add_option(
        "-y",
        "--ymax",
        dest="ymax",
        help="graphical y-axis limit for nucleotide misincorporation frequencies [%default]",
        type=float,
        default=0.3,
        action="store",
    )
    group2.add_option(
        "-m",
        "--readplot",
        dest="readplot",
        help="read length, in nucleotides, considered for plotting nucleotide misincorporations [%default]",
        type=int,
        default=25,
        action="store",
    )
    group2.add_option(
        "-b",
        "--refplot",
        dest="refplot",
        help="the number of reference nucleotides to consider for plotting base composition in the region located upstream "
        "and downstream of every read [%default]",
        type=int,
        default=10,
        action="store",
    )
    group2.add_option(
        "-t",
        "--title",
        dest="title",
        help="title used for plots [%default]",
        type="string",
        default="",
        action="store",
    )
    parser.add_option_group(group2)

    # Then the plethora of optional options for the statistical estimation ..
    group3 = OptionGroup(parser, "Options for the statistical estimation")
    group3.add_option(
        "",
        "--rand",
        dest="rand",
        help="Number of random starting points for the likelihood optimization  [%default]",
        type=int,
        default=30,
        action="store",
    )
    group3.add_option(
        "",
        "--burn",
        dest="burn",
        help="Number of burn-in iterations  [%default]",
        type=int,
        default=10000,
        action="store",
    )
    group3.add_option(
        "",
        "--adjust",
        dest="adjust",
        help="Number of adjust proposal variance parameters iterations  [%default]",
        type=int,
        default=10,
        action="store",
    )
    group3.add_option(
        "",
        "--iter",
        dest="iter",
        help="Number of final MCMC iterations  [%default]",
        type=int,
        default=50000,
        action="store",
    )
    group3.add_option(
        "",
        "--forward",
        dest="forward",
        help="Using only the 5' end of the seqs  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--reverse",
        dest="reverse",
        help="Using only the 3' end of the seqs  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--var-disp",
        dest="var_disp",
        help="Variable dispersion in the overhangs  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--jukes-cantor",
        dest="jukes_cantor",
        help="Use Jukes Cantor instead of HKY85  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--diff-hangs",
        dest="diff_hangs",
        help="The overhangs are different for 5' and 3'  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--fix-nicks",
        dest="fix_nicks",
        help="Fix the nick frequency vector (Only C.T from the 5' end and G.A from the 3' end)  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--use-raw-nick-freq",
        dest="use_raw_nick_freq",
        help="Use the raw nick frequency vector without smoothing  [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--single-stranded",
        dest="single_stranded",
        help="Single stranded protocol [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--theme-bw",
        dest="theme_bw",
        help="Use black and white theme in post. pred. plot [%default]",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "",
        "--seq-length",
        dest="seq_length",
        help="How long sequence to use from each side [%default]",
        type=int,
        default=12,
        action="store",
    )
    group3.add_option(
        "--stats-only",
        dest="stats_only",
        help="Run only statistical estimation from a valid result folder",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "--no-stats",
        help="Disabled statistical estimation, active by default",
        default=False,
        action="store_true",
    )
    group3.add_option(
        "--check-R-packages",
        help="Check if the R modules are working",
        default=False,
        action="store_true",
    )
    parser.add_option_group(group3)

    group4 = OptionGroup(parser, "Options for rescaling of BAM files")
    group4.add_option(
        "--rescale",
        dest="rescale",
        help="Rescale the quality scores in the BAM file using the output from the statistical estimation",
        default=False,
        action="store_true",
    )
    group4.add_option(
        "--rescale-only",
        dest="rescale_only",
        help="Run only rescaling from a valid result folder",
        default=False,
        action="store_true",
    )
    group4.add_option(
        "--rescale-out",
        dest="rescale_out",
        help="Write the rescaled BAM to this file",
        default=None,
        action="store",
    )
    group4.add_option(
        "--rescale-length-5p",
        dest="rescale_length_5p",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
        type=int,
        action="store",
    )
    group4.add_option(
        "--rescale-length-3p",
        dest="rescale_length_3p",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
        type=int,
        action="store",
    )
    parser.add_option_group(group4)

    # Parse the arguments
    (options, args) = parser.parse_args()

    # check python version
    if not check_py_version():
        return None

    # if the user wants to check the R packages then do that before the option parsing
    if options.check_R_packages:
        if check_R_lib():
            sys.exit(1)
        else:
            print("All R packages are present")
            sys.exit(0)

    # check general arguments
    if not (options.plot_only or options.stats_only) and not options.filename:
        parser.error("SAM/BAM file not given (-i)")
    if not (options.plot_only or options.ref):
        parser.error("Reference file not given (-r)")
    if not options.plot_only and not options.stats_only:
        if not file_exist(options.filename) or not file_exist(options.ref):
            return None
    if options.downsample is not None:
        if options.downsample <= 0:
            parser.error("-n/--downsample must be a positive value")
        elif options.downsample >= 1:
            options.downsample = int(options.downsample)

    if options.plot_only and not options.folder:
        parser.error("Folder not provided, required with --plot-only")
    if options.stats_only and not options.folder:
        parser.error("Folder not provided, required with --stats-only")
    if options.rescale_only and not options.folder:
        parser.error("Folder not provided, required with --rescale-only")
    if options.rescale_only and not options.filename:
        parser.error("Input bam not provided, required with --rescale-only")
    if options.rescale_only and not options.ref:
        parser.error("Reference not provided, required with --rescale-only")

    if options.verbose and options.quiet:
        parser.error("Cannot use verbose and quiet option at the same time")

    # check options
    if options.length < 0:
        parser.error("length (-l) must be a positive integrer")
    if options.around < 0:
        parser.error("around (-a) must be a positive integrer")
    if options.ymax <= 0 or options.ymax > 1:
        parser.error("ymax (-b) must be an real number between 0 and 1")
    if options.readplot < 0:
        parser.error("readplot (-m) must be a positive integrer")
    if options.refplot < 0:
        parser.error("refplot (-b) must be a positive integrer")
    if options.refplot > options.around and not options.plot_only:
        parser.error("refplot (-b) must be inferior to around (-a)")
    if options.readplot > options.length:
        parser.error("readplot (-m) must be inferior to length (-l)")
    if options.minqual < 0 or options.minqual > 41:
        parser.error(
            "minimal base quality, Phred score, must be within this range: 0 - 41"
        )

    # check statistic options
    if options.forward and options.reverse:
        parser.error(
            "Cannot use only forward end and only reverse end for the statistics"
        )

    # use filename as default for plot titles if not set
    if options.title == "" and options.filename:
        options.title = os.path.splitext(os.path.basename(options.filename))[0]
    # for --plot-only, use the folder name, without results_ as title
    if options.title == "" and not options.filename and options.folder:
        options.title = os.path.splitext(os.path.basename(options.folder))[0].replace(
            "results_", ""
        )

    # check folder
    if not options.folder and options.filename:
        options.folder = (
            "results_" + os.path.splitext(os.path.basename(options.filename))[0]
        )

    # check destination for rescaled bam
    if not options.rescale_out and (options.rescale or options.rescale_only):
        # if there are multiple bam files to rescale then pick first one as
        # the name of the rescaled file
        if isinstance(options.filename, list):
            basename = os.path.basename(options.filename[0])
        else:
            basename = os.path.basename(options.filename)
        with_ext = os.path.splitext(basename)[0] + ".rescaled.bam"
        options.rescale_out = os.path.join(options.folder, with_ext)

    if os.path.isdir(options.folder):
        if not options.quiet and not options.plot_only:
            print(("Warning, %s already exists" % options.folder))
        if options.plot_only:
            if not file_exist(options.folder + "/dnacomp.txt") or not file_exist(
                options.folder + "/misincorporation.txt"
            ):
                parser.error("folder %s is not a valid result folder" % options.folder)
    else:
        os.makedirs(options.folder, mode=0o750)
        if options.plot_only or options.stats_only or options.rescale_only:
            sys.stderr.write(
                "Error, %s does not exist while plot/stats/rescale only was used\n"
                % options.folder
            )
            return None

    if options.rescale_length_3p is None:
        options.rescale_length_3p = options.seq_length
    elif not (0 <= options.rescale_length_3p <= options.seq_length):
        parser.error(
            "--rescale-length-3p must be less than or equal to "
            "--seq-length and greater than zero"
        )

    if options.rescale_length_5p is None:
        options.rescale_length_5p = options.seq_length
    elif not (0 <= options.rescale_length_5p <= options.seq_length):
        parser.error(
            "--rescale-length-5p must be less than or equal to "
            "--seq-length and greater than zero"
        )

    # check if the Rscript executable is present on the system
    if not whereis("Rscript"):
        print("Warning, Rscript is not in your PATH, plotting is disabled")
        options.no_r = True

    # check the nick frequencies options
    if (options.use_raw_nick_freq + options.fix_nicks + options.single_stranded) > 1:
        parser.error(
            "The options --use-raw-nick-freq, --fix-nicks and --single-stranded are mutually exclusive."
        )

    if check_R_lib():
        # check for R libraries
        print("The Bayesian estimation has been disabled\n")
        options.no_stats = True
        if options.stats_only:
            sys.exit("Cannot use --stats-only with missing R libraries")
        if options.rescale:
            sys.exit("Cannot use --rescale with missing R libraries")
        if options.rescale_only:
            sys.exit("Cannot use --rescale-only with missing R libraries")

    return options
