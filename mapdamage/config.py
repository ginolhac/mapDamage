import argparse
import logging
import shutil
import sys

from argparse import ArgumentError
from pathlib import Path

from mapdamage.version import __version__
from mapdamage.rscript import check_r_libraries


def file_exist(filename):
    return filename == Path("-") or filename.is_file()


class NumericParser:
    """Type checking for numerical arguments within some range"""

    def __init__(self, cls, min_value=float("-inf"), max_value=float("inf")):
        self._cls = cls
        self._min = min_value
        self._max = max_value

    def __call__(self, value):
        value = self._cls(value)
        if value < self._min:
            raise argparse.ArgumentTypeError(
                "must be greater than or equal to %s" % (self._min,)
            )
        elif value > self._max:
            raise argparse.ArgumentTypeError(
                "must be less than or equal to %s" % (self._max,)
            )

        return value

    def __repr__(self):
        return "numeric"


UnsignedInt = NumericParser(int, min_value=0)
PositiveInt = NumericParser(int, min_value=1)


class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Modified ArgumentDefaultsHelpFormatter that excludes several constants (True,
    False, None) and uses a custom presentation of the default value.
    """

    def _get_help_string(self, action):
        # The following values look silly as part of a help string
        if action.default is None or isinstance(action.default, bool):
            return action.help

        # The subclass does not allow modification to the defaults string, so instead
        # we access the logic by simply checking if the result was modified.
        if super()._get_help_string(action) == action.help:
            return action.help

        return action.help + " [%(default)s]"


class ThrowingArgumentParser(argparse.ArgumentParser):
    def exit(self, status=0, message=None):
        if status:
            raise ArgumentError(None, message.strip() if message else None)

        sys.exit(status)

    def error(self, message):
        # Workaround for hard-coded formatting of some errors
        error = sys.exc_info()[1]
        if isinstance(error, ArgumentError):
            raise error

        self.exit(2, message)


def _build_parser():
    parser = ThrowingArgumentParser(
        prog="mapDamage",
        usage="%(prog)s [options] -i alignment.bam -r reference.fasta",
        epilog="Please report bugs on GitHub: "
        "https://github.com/ginolhac/mapDamage/issues/new",
        formatter_class=CustomHelpFormatter,
    )

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    group = parser.add_argument_group("Input and output")
    group.add_argument(
        "-i",
        "--input",
        dest="filename",
        type=Path,
        metavar="SAM/BAM",
        help="SAM/BAM file; use '-' to read from stdin",
    )
    group.add_argument(
        "-r",
        "--reference",
        dest="ref",
        type=Path,
        metavar="FASTA",
        help="Reference genome in FASTA format",
    )
    group.add_argument(
        "-d",
        "--folder",
        help="Path of folder used to store results [results_FILENAME]",
        type=Path,
    )
    group.add_argument(
        "-n",
        "--downsample",
        type=float,
        metavar="X",
        help="Downsample to a randomly selected fraction of the reads (if 0 < "
        "DOWNSAMPLE < 1), or a fixed number of randomly selected reads (if DOWNSAMPLE "
        ">= 1). By default, no downsampling is performed.",
    )
    group.add_argument(
        "--downsample-seed",
        type=int,
        metavar="X",
        help="Seed value to use for downsampling. See documentation for py module "
        "'random' for default behavior.",
    )

    group = parser.add_argument_group("General options")
    group.add_argument(
        "--merge-libraries",
        action="store_true",
        help="Treat BAM as containing only a single library",
    )
    group.add_argument(
        "--merge-reference-sequences", help=argparse.SUPPRESS, action="store_true"
    )
    group.add_argument(
        "-l",
        "--length",
        type=PositiveInt,
        default=70,
        metavar="N",
        help="Number of nucleotides to process from the 5' and 3' of alignments",
    )
    group.add_argument(
        "-a",
        "--around",
        type=UnsignedInt,
        default=10,
        metavar="N",
        help="Number of nucleotides to process before/after alignments",
    )
    group.add_argument(
        "-Q",
        "--min-basequal",
        dest="minqual",
        type=NumericParser(int, 0, 93),
        default=0,
        metavar="PHRED",
        help="minimun base quality Phred score considered, Phred-33 assumed",
    )
    group.add_argument(
        "--plot-only",
        help="Run only plotting from a valid result folder",
        action="store_true",
    )
    group.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        metavar="LEVEL",
        type=str.upper,
        help="Logging verbosity level; one of DEBUG, INFO, WARNING, and ERROR",
    )
    group.add_argument(
        "--no-plot",
        dest="no_r",
        help=argparse.SUPPRESS,
        action="store_true",
    )

    # options for plotting damage patterns
    group = parser.add_argument_group("Options for plots")
    group.add_argument(
        "-y",
        "--ymax",
        type=float,
        default=0.3,
        metavar="Y",
        help="Upper limit for y-axis in nucleotide misincorporation plots",
    )
    group.add_argument(
        "-m",
        "--readplot",
        type=PositiveInt,
        default=25,
        metavar="N",
        help="Number of bases plotted from the 5' and 3' termini ead length, in "
        "nucleotides, considered for plotting nucleotide misincorporations",
    )
    group.add_argument(
        "-b",
        "--refplot",
        type=PositiveInt,
        default=10,
        metavar="N",
        help="the number of reference nucleotides to consider for ploting base "
        "composition in the region located upstream and downstream of every read",
    )
    group.add_argument("-t", "--title", help="title used for plots")

    # Then the plethora of optional options for the statistical estimation ..
    group = parser.add_argument_group("Options for the statistical estimation")
    group.add_argument(
        "--rand",
        type=PositiveInt,
        default=30,
        metavar="X",
        help="Number of random starting points for the likelihood optimization",
    )
    group.add_argument(
        "--burn",
        type=PositiveInt,
        default=10000,
        metavar="N",
        help="Number of burn-in iterations",
    )
    group.add_argument(
        "--adjust",
        type=int,
        default=10,
        metavar="N",
        help="Number of adjust proposal variance parameter iterations",
    )
    group.add_argument(
        "--iter",
        type=PositiveInt,
        default=50000,
        metavar="N",
        help="Number of final MCMC iterations",
    )
    group.add_argument(
        "--termini",
        default="both",
        choices=("5p", "3p", "both"),
        help="Use either mismatches at 5p (forward), 3p (reverse), or both termini for "
        "damage models",
    )
    group.add_argument(
        "--forward",
        action="store_const",
        const="5p",
        dest="termini",
        help=argparse.SUPPRESS,
    )
    group.add_argument(
        "--reverse",
        action="store_const",
        const="3p",
        dest="termini",
        help=argparse.SUPPRESS,
    )
    group.add_argument(
        "--var-disp",
        action="store_true",
        help="Variable dispersion in the overhangs",
    )
    group.add_argument(
        "--jukes-cantor",
        action="store_true",
        help="Use Jukes Cantor instead of HKY85",
    )
    group.add_argument(
        "--diff-hangs",
        action="store_true",
        help="The overhangs are different for 5' and 3'",
    )
    group.add_argument(
        "--fix-nicks",
        action="store_true",
        help="Fix the nick frequency vector: Only C>T from the 5' end and G>A from the "
        "3' end)",
    )
    group.add_argument(
        "--use-raw-nick-freq",
        action="store_true",
        help="Use the raw nick frequency vector without smoothing",
    )
    group.add_argument(
        "--single-stranded",
        action="store_true",
        help="Single stranded protocol",
    )
    group.add_argument(
        "--theme-bw",
        action="store_true",
        help="Use black and white theme in post. pred. plot",
    )
    group.add_argument(
        "--seq-length",
        type=int,
        default=12,
        metavar="N",
        help="How long sequence to use from each side",
    )
    group.add_argument(
        "--stats-only",
        action="store_true",
        help="Run only statistical estimation. Requires a valid --folder",
    )
    group.add_argument(
        "--no-stats",
        action="store_true",
        help="Disabled statistical estimation, active by default",
    )
    group.add_argument(
        "--check-R-packages",
        action="store_true",
        help="Check if required R modules are available; terminates after checking",
    )

    group = parser.add_argument_group("Options for rescaling of BAM files")
    group.add_argument(
        "--rescale",
        action="store_true",
        help="Rescale the quality scores in the BAM file using the output from the "
        "statistical estimation",
    )
    group.add_argument(
        "--rescale-only",
        action="store_true",
        help="Run only rescaling. Requires a valid --folder and --input file",
    )
    group.add_argument(
        "--rescale-out",
        metavar="BAM",
        help="Write the rescaled BAM to this file. By default, the rescaled BAM saved "
        "in --folder, using the same name as the --input file but with a .rescaled.bam "
        "extension.",
    )
    group.add_argument(
        "--rescale-length-5p",
        type=int,
        metavar="N",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
    )
    group.add_argument(
        "--rescale-length-3p",
        type=int,
        metavar="N",
        help="How many bases to rescale at the 5' termini; defaults to --seq-length.",
    )

    return parser


def parse_args(argv):
    parser = _build_parser()
    options = parser.parse_args(argv)
    logger = logging.getLogger(__name__)

    # Set logging levels for root logger and any additional handlers
    logging.getLogger().setLevel(options.log_level)
    for handler in logging.getLogger().handlers:
        handler.setLevel(options.log_level)

    # check if the Rscript executable is present on the system
    if not shutil.which("Rscript"):
        logger.warning("Rscript is not in your PATH, plotting is disabled")
        options.no_r = True

    # if the user wants to check the R packages then do that before the option parsing
    if options.check_R_packages:
        if options.no_r:
            logger.error("Cannot check for R packages without Rscript")
            parser.exit(1)
        elif not check_r_libraries():
            parser.exit(1)
        else:
            logger.info("All R packages are present")
            parser.exit()

    # check general arguments
    if not (options.plot_only or options.stats_only) and not options.filename:
        parser.error("--input SAM/BAM file not specified")
    if not (options.plot_only or options.ref):
        parser.error("--reference FASTA file not specified")
    if not options.plot_only and not options.stats_only:
        if not file_exist(options.filename) or not file_exist(options.ref):
            logger.error("%s is not a valid file", options.filename)
    if options.downsample is not None:
        if options.downsample <= 0:
            parser.error("-n/--downsample must be a positive value")
        elif options.downsample >= 1:
            options.downsample = int(options.downsample)

    if options.plot_only and not options.folder:
        parser.error("--folder required when using --plot-only")
    if options.stats_only and not options.folder:
        parser.error("--folder required when using --stats-only")
    if options.rescale_only:
        if not options.folder:
            parser.error("--folder required when using --rescale-only")
        elif not options.filename:
            parser.error("--input required when using --rescale-only")
        elif not options.ref:
            parser.error("--reference required when using --rescale-only")

    # check options
    if options.ymax <= 0 or options.ymax > 1:
        parser.error("--ymax (-b) must be an real number beetween 0 and 1")
    if options.refplot > options.around and not options.plot_only:
        parser.error("--refplot (-b) must be less than --around (-a)")
    if options.readplot > options.length:
        parser.error("--readplot (-m) must be less than --length (-l)")

    # use filename as default for plot titles if not set
    if options.title is None:
        if options.filename:
            options.title = options.filename.stem
        # for --plot-only, use the folder name, without results_ as title
        elif options.folder:
            options.title = options.folder.stem.replace("results_", "")
        else:
            options.title = ""

    # check folder
    if not options.folder and options.filename:
        options.folder = Path(options.filename.stem + ".mapDamage")

    # check destination for rescaled bam
    if not options.rescale_out and (options.rescale or options.rescale_only):
        options.rescale_out = options.folder / (options.filename.stem + ".rescaled.bam")

    if options.folder.is_dir():
        if not options.plot_only:
            logger.warning(
                "Folder '%s' already exists; content may be overwritten", options.folder
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
            "The options --use-raw-nick-freq, --fix-nicks and --single-stranded are "
            "mutually exclusive."
        )

    if options.no_r or not check_r_libraries():
        # check for R libraries
        logger.warning("The Bayesian estimation has been disabled")
        options.no_stats = True
        if options.stats_only:
            parser.error("Cannot use --stats-only with missing R libraries")
        if options.rescale:
            parser.error("Cannot use --rescale with missing R libraries")
        if options.rescale_only:
            parser.error("Cannot use --rescale-only with missing R libraries")

    return options
