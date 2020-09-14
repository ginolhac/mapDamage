import logging
import os
import subprocess
import time

from pathlib import Path
from pkg_resources import resource_filename

from mapdamage.version import __version__


def misincorporation_plot(options):
    folder = options.folder.absolute()
    fmut = folder / "misincorporation.txt"
    fcomp = folder / "dnacomp.txt"
    output = folder / "Fragmisincorporation_plot.pdf"

    logger = logging.getLogger(__name__)
    logger.info("Saving misincorporation plot to '%s'", output)

    return _rscript_call(
        Path("mapDamage.r"),
        COMP=fcomp,
        PDFOUT=output,
        AROUND=options.refplot,
        MISINCORP=fmut,
        LENGTH=options.readplot,
        YMAX=options.ymax,
        FOLDER=folder,
        TITLE=options.title,
        VERSION=__version__,
    )


def length_distribution_plot(options):
    """optional length distribution and cumulative C>T mutations plots, per strand"""
    folder = options.folder.absolute()
    fmut = folder / "misincorporation.txt"
    flength = folder / "lgdistribution.txt"
    output = folder / "Length_plot.pdf"

    logger = logging.getLogger(__name__)
    logger.info("Saving length distribution plot to '%s'", output)

    return _rscript_call(
        Path("lengths.r"),
        LGDIST=flength,
        PDFOUT=output,
        MISINCORP=fmut,
        TITLE=options.title,
        VERSION=__version__,
    )


def check_r_libraries():
    """Checks if the necessary R libraries are here, signal otherwise"""
    logger = logging.getLogger(__name__)
    missing_libries = False

    for library in ["ggplot2", "gam", "Rcpp", "RcppGSL"]:
        command = ["Rscript", "-e", "library(%s)" % (library,)]

        if not _log_call(command, quiet=True):
            logger.error("Required R library is missing: %r", library)
            missing_libries = True

    return not missing_libries


def perform_bayesian_estimates(options):
    """Runs the Bayesian estimation program"""
    logger = logging.getLogger(__name__)
    logger.info("Performing Bayesian estimates")
    folder = options.folder.absolute()

    # Disable compile time warnings for Rcpp code; these are outside of our control
    env = dict(os.environ)
    env["PKG_CPPFLAGS"] = env.get("PKG_CPPFLAGS", "") + " -w"

    return _rscript_call(
        Path("stats") / "runGeneral.r",
        GRID_ITER=options.rand,
        BURN_IN=options.burn,
        ADJUST_ITER=options.adjust,
        ITERATIONS=options.iter,
        TERMINI=options.termini,
        FIX_DISP=not options.var_disp,
        SAME_OVERHANGS=not options.diff_hangs,
        FIX_NU=bool(options.fix_nicks),
        DS_PROTOCOL=not options.single_stranded,
        SUB_LENGTH=options.seq_length,
        PATH_TO_DAT=str(folder) + "/",
        VERBOSE=bool(options.log_level == "DEBUG"),
        QUIET=bool(options.log_level not in ("DEBUG", "INFO")),
        JUKES_CANTOR=bool(options.jukes_cantor),
        USE_RAW_NICK_FREQ=bool(options.use_raw_nick_freq),
        USE_BW_THEME=bool(options.theme_bw),
        env=env,
    )


def _rscript_call(filepath, env=None, **kwargs):
    cwd = Path(resource_filename("mapdamage", "r")) / filepath.parent
    command = ["Rscript", filepath.name]
    for item in sorted(kwargs.items()):
        command.append("%s=%s" % item)

    return _log_call(command, cwd=cwd, env=env)


def _log_call(command, quiet=False, cwd=None, env=None):
    command = [str(value) for value in command]

    logger = logging.getLogger(__name__)
    logger.debug("Running command %r", " ".join(command))
    loglevel = logging.DEBUG if quiet else logging.INFO

    start = time.time()

    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        cwd=cwd,
        env=env,
    )

    try:
        for line in proc.stdout:
            logger.log(loglevel, "%s", line.decode("utf-8", errors="replace").rstrip())

        returncode = proc.wait()
    except:
        proc.terminate()
        proc.wait()
        raise

    logger.debug("Call completed in %.2fs with rc %i", time.time() - start, returncode)

    if returncode and not quiet:
        logger.error("Command returned error %i: %r", returncode, " ".join(command))

    return not returncode
