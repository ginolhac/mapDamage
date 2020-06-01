import logging
import os
import subprocess
import time

from pkg_resources import resource_filename

from mapdamage.version import __version__


def construct_path(name, folder="Rscripts"):
    """Construct a path to the mapdamage package data given the name
    """
    return resource_filename("mapdamage", os.path.join(folder, name))


def misincorporation_plot(opt):
    fmut = opt.folder / "misincorporation.txt"
    fcomp = opt.folder / "dnacomp.txt"
    output = opt.folder / "Fragmisincorporation_plot.pdf"

    logger = logging.getLogger(__name__)
    logger.info("Saving misincorporation plot to %r", output)

    return _log_call(
        [
            "Rscript",
            construct_path("mapDamage.R"),
            fcomp,
            output,
            opt.refplot,
            fmut,
            opt.readplot,
            opt.ymax,
            opt.folder,
            opt.title,
            __version__,
        ]
    )


def length_distribution_plot(opt):
    """optional length distribution and cumulative C>T mutations plots, per strand
    """
    fmut = opt.folder / "misincorporation.txt"
    flength = opt.folder / "lgdistribution.txt"
    output = opt.folder / "Length_plot.pdf"

    logger = logging.getLogger(__name__)
    logger.info("Saving length distribution plot to %r", output)

    return _log_call(
        [
            "Rscript",
            construct_path("lengths.R"),
            flength,
            output,
            fmut,
            opt.length,
            opt.title,
            __version__,
            int(opt.log_level not in ("DEBUG", "INFO")),
        ]
    )


def check_r_libraries():
    """Checks if the necessary R libraries are here, signal otherwise
    """
    logger = logging.getLogger(__name__)
    script = construct_path("stats/checkLibraries.R")
    missing_libries = False
    with open(os.devnull, "w") as null:
        for library in ["ggplot2", "gam", "Rcpp", "RcppGSL"]:
            command = ["Rscript", script, library]

            if not _log_call(command, log_failures=False, stdout=null, stderr=null):
                logger.error("Required R library is missing: %r", library)
                missing_libries = True

    return not missing_libries


def perform_bayesian_estimates(opt):
    """Runs the Bayesian estimation program, using the options opt
    """
    logger = logging.getLogger(__name__)
    logger.info("Performing Bayesian estimates")

    return _log_call(
        [
            "Rscript",
            construct_path("stats/runGeneral.R"),
            opt.rand,
            opt.burn,
            opt.adjust,
            opt.iter,
            int(opt.forward),
            int(opt.reverse),
            int(not opt.var_disp),
            int(not opt.diff_hangs),
            0,
            int(opt.fix_nicks),
            int(not opt.single_stranded),
            opt.seq_length,
            str(opt.folder) + "/",
            construct_path("stats/"),
            opt.folder / "Stats_out",
            int(opt.log_level == "DEBUG"),
            int(opt.log_level not in ("DEBUG", "INFO")),
            int(opt.jukes_cantor),
            opt.folder / "acgt_ratio.csv",
            int(opt.use_raw_nick_freq),
            int(opt.theme_bw),
        ]
    )


def _log_call(command, log_failures=True, **kwargs):
    command = [str(value) for value in command]

    logger = logging.getLogger(__name__)
    logger.debug("Running command %r", " ".join(command))

    start = time.time()
    returncode = subprocess.call(command, **kwargs)
    logger.debug("Call completed in %.2fs with rc %i", time.time() - start, returncode)

    if returncode and log_failures:
        logger.error("Command returned error %i: %r", returncode, " ".join(command))

    return not returncode
