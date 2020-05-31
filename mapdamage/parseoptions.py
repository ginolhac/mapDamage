import logging
import shutil
import sys

from argparse import ArgumentParser, SUPPRESS
from pathlib import Path

from mapdamage.version import __version__
from mapdamage.rscript import check_r_libraries


def file_exist(filename):
    return filename == Path("-") or filename.is_file()


def _build_parser():
    parser = ArgumentParser(
        prog="mapDamage",
        usage="%(prog)s [options] -i BAMfile -r reference.fasta",
        epilog="Please report bugs on GitHub: https://github.com/ginolhac/mapDamage/issues/new",
    )

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.add_argument_group("Input files")
    args.add_argument(
        "-i",
        "--input",
        help="SAM/BAM file, must contain a valid header, use '-' for reading a BAM from stdin",
        dest="filename",
        type=Path,
    )
    args.add_argument(
        "-r",
        "--reference",
        help="Reference file in FASTA format",
        dest="ref",
        type=Path,
    )

    parser.add_argument_group(args)
    group = parser.add_argument_group("General options")
    group.add_argument(
        "-n",
        "--downsample",
        help="Downsample to a randomly selected fraction of the reads (if 0 < DOWNSAMPLE < 1), or "
        "a fixed number of randomly selected reads (if DOWNSAMPLE >= 1). By default, no downsampling is performed.",
        type=float,
    )
    group.add_argument(
        "--downsample-seed",
        help="Seed value to use for downsampling. See documentation for py module 'random' for default behavior.",
        type=int,
    )
    group.add_argument(
        "--merge-reference-sequences",
        help=SUPPRESS,
        default=False,
        action="store_true",
    )
    group.add_argument(
        "-l",
        "--length",
        help="read length, in nucleotides to consider [%(default)s]",
        type=int,
        default=70,
    )
    group.add_argument(
        "-a",
        "--around",
        help="nucleotides to retrieve before/after reads [%(default)s]",
        type=int,
        default=10,
    )
    group.add_argument(
        "-Q",
        "--min-basequal",
        dest="minqual",
        help="minimun base quality Phred score considered, Phred-33 assumed [%(default)s]",
        type=int,
        default=0,
    )
    group.add_argument(
        "-d",
        "--folder",
        help="folder name to store results [results_FILENAME]",
        type=Path,
    )
    group.add_argument(
        "--plot-only",
        help="Run only plotting from a valid result folder",
        default=False,
        action="store_true",
    )
    group.add_argument(
        "--log-level",
        help="Logging verbosity level; one of DEBUG, INFO, WARNING, and ERROR [%(default)s]",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    group.add_argument(
        "--no-plot", dest="no_r", help=SUPPRESS, default=False, action="store_true"
    )
    parser.add_argument_group(group)

    # options for plotting damage patterns
    group2 = parser.add_argument_group("Options for graphics")
    group2.add_argument(
        "-y",
        "--ymax",
        help="graphical y-axis limit for nucleotide misincorporation frequencies [%(default)s]",
        type=float,
        default=0.3,
    )
    group2.add_argument(
        "-m",
        "--readplot",
        help="read length, in nucleotides, considered for plotting nucleotide misincorporations [%(default)s]",
        type=int,
        default=25,
    )
    group2.add_argument(
        "-b",
        "--refplot",
        help="the number of reference nucleotides to consider for ploting base composition in the region located upstream "
        "and downstream of every read [%(default)s]",
        type=int,
        default=10,
    )
    group2.add_argument(
        "-t", "--title", help="title used for plots [%(default)s]", default="",
    )
    parser.add_argument_group(group2)

    # Then the plethora of optional options for the statistical estimation ..
    group3 = parser.add_argument_group("Options for the statistical estimation")
    group3.add_argument(
        "--rand",
        help="Number of random starting points for the likelihood optimization [%(default)s]",
        type=int,
        default=30,
    )
    group3.add_argument(
        "--burn",
        help="Number of burnin iterations [%(default)s]",
        type=int,
        default=10000,
    )
    group3.add_argument(
        "--adjust",
        help="Number of adjust proposal variance parameters iterations [%(default)s]",
        type=int,
        default=10,
    )
    group3.add_argument(
        "--iter",
        help="Number of final MCMC iterations [%(default)s]",
        type=int,
        default=50000,
    )
    group3.add_argument(
        "--forward",
        help="Using only the 5' end of the seqs [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--reverse",
        help="Using only the 3' end of the seqs [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--var-disp",
        help="Variable dispersion in the overhangs [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--jukes-cantor",
        help="Use Jukes Cantor instead of HKY85 [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--diff-hangs",
        help="The overhangs are different for 5' and 3' [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--fix-nicks",
        help="Fix the nick frequency vector (Only C.T from the 5' end and G.A from the 3' end) [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--use-raw-nick-freq",
        help="Use the raw nick frequency vector without smoothing [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--single-stranded",
        help="Single stranded protocol [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--theme-bw",
        help="Use black and white theme in post. pred. plot [%(default)s]",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--seq-length",
        help="How long sequence to use from each side [%(default)s]",
        type=int,
        default=12,
    )
    group3.add_argument(
        "--stats-only",
        help="Run only statistical estimation from a valid result folder",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--no-stats",
        help="Disabled statistical estimation, active by default",
        default=False,
        action="store_true",
    )
    group3.add_argument(
        "--check-R-packages",
        help="Check if the R modules are working",
        default=False,
        action="store_true",
    )
    parser.add_argument_group(group3)

    group4 = parser.add_argument_group("Options for rescaling of BAM files")
    group4.add_argument(
        "--rescale",
        help="Rescale the quality scores in the BAM file using the output from the statistical estimation",
        default=False,
        action="store_true",
    )
    group4.add_argument(
        "--rescale-only",
        help="Run only rescaling from a valid result folder",
        default=False,
        action="store_true",
    )
    group4.add_argument(
        "--rescale-out", help="Write the rescaled BAM to this file",
    )
    group4.add_argument(
        "--rescale-length-5p",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
        type=int,
    )
    group4.add_argument(
        "--rescale-length-3p",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
        type=int,
    )
    parser.add_argument_group(group4)

    return parser


def options(argv):
    parser = _build_parser()
    options = parser.parse_args(argv)
    logger = logging.getLogger(__name__)

    # check if the Rscript executable is present on the system
    if not shutil.which("Rscript"):
        logger.warning("Rscript is not in your PATH, plotting is disabled")
        options.no_r = True

    # if the user wants to check the R packages then do that before the option parsing
    if options.check_R_packages:
        if options.no_r:
            logger.error("Cannot check for R packages without Rscript")
            sys.exit(1)
        elif not check_r_libraries():
            sys.exit(1)
        else:
            logger.info("All R packages are present")
            sys.exit(0)

    # check general arguments
    if not (options.plot_only or options.stats_only) and not options.filename:
        parser.error("SAM/BAM file not given (-i)")
    if not (options.plot_only or options.ref):
        parser.error("Reference file not given (-r)")
    if not options.plot_only and not options.stats_only:
        if not file_exist(options.filename) or not file_exist(options.ref):
            logger.error("%s is not a valid file", options.filename)
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

    # check options
    if options.length < 0:
        parser.error("length (-l) must be a positive integrer")
    if options.around < 0:
        parser.error("around (-a) must be a positive integrer")
    if options.ymax <= 0 or options.ymax > 1:
        parser.error("ymax (-b) must be an real number beetween 0 and 1")
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
        options.title = options.filename.stem
    # for --plot-only, use the folder name, without results_ as title
    if options.title == "" and not options.filename and options.folder:
        options.title = options.folder.stem.replace("results_", "")

    # check folder
    if not options.folder and options.filename:
        options.folder = "results_" + options.filename.stem

    # check destination for rescaled bam
    if not options.rescale_out and (options.rescale or options.rescale_only):
        options.rescale_out = options.folder / (options.filename.stem + ".rescaled.bam")

    if options.folder.is_dir():
        if not options.plot_only:
            logger.warning(
                "Folder %r already exists; content may be overwritten", options.folder
            )
        if options.plot_only:
            if not (
                (options.folder / "dnacomp.txt").is_file()
                and (options.folder / "misincorporation.txt").is_file()
            ):
                parser.error("'%s' is not a valid result folder" % options.folder)
    else:
        options.folder.mkdir(parents=True, exist_ok=True, mode=0o750)
        if options.plot_only or options.stats_only or options.rescale_only:
            logger.error(
                "Folder %s does not exist while plot/stats/rescale only was used",
                options.folder,
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

    # check the nick frequencies options
    if (options.use_raw_nick_freq + options.fix_nicks + options.single_stranded) > 1:
        parser.error(
            "The options --use-raw-nick-freq, --fix-nicks and --single-stranded are mutually exclusive."
        )

    if options.no_r or not check_r_libraries():
        # check for R libraries
        logger.warning("The Bayesian estimation has been disabled")
        options.no_stats = True
        if options.stats_only:
            sys.exit("Cannot use --stats-only with missing R libraries")
        if options.rescale:
            sys.exit("Cannot use --rescale with missing R libraries")
        if options.rescale_only:
            sys.exit("Cannot use --rescale-only with missing R libraries")

    return options
