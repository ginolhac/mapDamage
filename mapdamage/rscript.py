import logging
import os
import subprocess
import time

from subprocess import CalledProcessError, check_call

import mapdamage

from mapdamage.version import __version__


def construct_path(name, folder="Rscripts"):
    """Construct a path to the mapdamage package data given the name
    """
    return os.path.join(mapdamage.__path__[0], folder, name)


def plot(opt):
    """Calls the 'Rscript' to draw the plots.
    """
    fmut = opt.folder + "/" + "misincorporation.txt"
    fcomp = opt.folder + "/" + "dnacomp.txt"
    title = opt.folder + "/" + "Fragmisincorporation_plot.pdf"

    script = construct_path("mapDamage.R")
    call = [
        "Rscript",
        script,
        fcomp,
        title,
        opt.refplot,
        fmut,
        opt.readplot,
        opt.ymax,
        opt.folder,
        opt.title,
        __version__,
    ]
    code = subprocess.call(list(map(str, call)))

    logger = logging.getLogger(__name__)
    if code == 0:
        logger.info("pdf %s generated", title)
        return 0
    else:
        logger.error("plotting with R failed")
        return 1


def opt_plots(opt):
    """optional length distribution and cumulative C>T mutations plots, per strand
    """
    fmut = opt.folder + "/" + "misincorporation.txt"
    flength = opt.folder + "/" + "lgdistribution.txt"
    output = opt.folder + "/" + "Length_plot.pdf"

    script = construct_path("lengths.R")
    call = [
        "Rscript",
        script,
        flength,
        output,
        fmut,
        opt.length,
        opt.title,
        __version__,
        int(opt.log_level not in ("DEBUG", "INFO")),
    ]
    code = subprocess.call(list(map(str, call)))
    if code:
        logger = logging.getLogger(__name__)
        logger.error("plotting with R failed with return-code %i", code)

    return code


def check_R_one_lib(name):
    """Checks if a necessary R library is here
    """
    with open(os.devnull, "w") as devnull:
        rpa = construct_path("stats/checkLibraries.R")
        return subprocess.call(
            ["Rscript", rpa, "--args", name], stdout=devnull, stderr=devnull
        )


def check_R_lib():
    """Checks if the necessary R libraries are here, signal otherwise
    """
    libs = ["inline", "ggplot2", "gam", "Rcpp", "RcppGSL"]
    missing_lib = []
    for lib in libs:
        # Check the libraries
        if check_R_one_lib(lib):
            # found a missing library
            missing_lib.append(lib)

    logger = logging.getLogger(__name__)
    if len(missing_lib) > 1:
        # Grammar Nazi has arrived
        last_ele = missing_lib.pop()
        logger.error(
            "Missing the following R libraries '"
            + "', '".join(missing_lib)
            + "' and '"
            + last_ele
            + "'"
        )
        return 1
    elif len(missing_lib) == 1:
        logger.error("Missing the following R library %s", missing_lib[0])
        return 1
    else:
        # No missing libraries
        return 0


def run_stats(opt):
    """Runs the Bayesian estimation program, using the options opt
    """
    arg = [
        "Rscript",
        construct_path("stats/runGeneral.R"),
        "--args",
        opt.rand,
        opt.burn,
        opt.adjust,
        opt.iter,
        bool_to_int(opt.forward),
        bool_to_int(opt.reverse),
        bool_to_int(not opt.var_disp),
        bool_to_int(not opt.diff_hangs),
        0,
        bool_to_int(opt.fix_nicks),
        bool_to_int(not opt.single_stranded),
        opt.seq_length,
        opt.folder + "/",
        construct_path("stats/"),
        opt.folder + "/Stats_out",
        1 if opt.log_level == "DEBUG" else 0,
        1 if opt.log_level not in ("DEBUG", "INFO") else 0,
        int(opt.jukes_cantor),
        opt.folder + "/acgt_ratio.csv",
        int(opt.use_raw_nick_freq),
        int(opt.theme_bw),
    ]
    arg = [str(i) for i in arg]

    logger = logging.getLogger(__name__)
    logger.info("Performing Bayesian estimates")
    logger.debug("Call: %s", " ".join(arg))
    start_time = time.time()

    try:
        check_call(arg)
    except CalledProcessError as e:
        logger.error("The Bayesian statistics program failed to finish")
        raise e

    logger.debug("Bayesian estimates completed in %f seconds", time.time() - start_time)
    return 0


def bool_to_int(boolean):
    """ return 0 for False and 1 for True """
    result = 1 if boolean else 0
    return result
