# =============================================================================
# 00_setup.R  --  locate the Python engine and the conda environment
# =============================================================================
#
# Every other script in R/ assumes this file has been sourced first. It does
# two things:
#
#   1. remembers where the repository lives, so that scripts can be run from
#      any working directory;
#   2. connects reticulate to the conda environment that holds msprime.
#
# The conda environment is *not* hard-coded. It is taken, in order, from
#   - the argument you pass to init_python()
#   - the environment variable GLACIALSIM_CONDA_ENV
#   - the default name "msprime-env"
#
# So on a new machine you either do
#     Sys.setenv(GLACIALSIM_CONDA_ENV = "/home/me/miniconda3/envs/msprime-env")
# once per session (or in ~/.Renviron), or call
#     init_python("/home/me/miniconda3/envs/msprime-env")
# =============================================================================

suppressPackageStartupMessages(library(reticulate))

## ---- repository root --------------------------------------------------------
# Set GLACIALSIM_HOME if you source these files from somewhere unusual;
# otherwise the working directory is assumed to be the repository root.
glacialsim_home <- function() {
  h <- Sys.getenv("GLACIALSIM_HOME", unset = "")
  if (nzchar(h)) normalizePath(h, mustWork = TRUE) else normalizePath(getwd())
}

glacialsim_python_dir <- function() file.path(glacialsim_home(), "python")


## ---- Python / msprime -------------------------------------------------------
# Connects reticulate to the conda environment and imports python/glacial_sim.py
# as a namespaced module object. Returns that module, and also stores it in the
# global variable `gs` so the wrappers in 01_simulate.R can find it.
#
# conda_env : name of a conda environment ("msprime-env") or a full path to it.
# quiet     : suppress the version report.
init_python <- function(conda_env = NULL, quiet = FALSE) {

  if (is.null(conda_env))
    conda_env <- Sys.getenv("GLACIALSIM_CONDA_ENV", unset = "msprime-env")

  reticulate::use_condaenv(conda_env, required = TRUE)

  gs <<- reticulate::import_from_path("glacial_sim",
                                      path = glacialsim_python_dir(),
                                      convert = TRUE)

  # The engine hands numpy arrays back to R. If reticulate cannot convert them
  # this fails later, deep inside a wrapper, with an unhelpful message about
  # coercing an environment - so check it here instead. The usual cause is an
  # old reticulate paired with numpy 2.x: reticulate gained numpy 2 support in
  # 1.38, and msprime 1.4 onwards requires numpy >= 2.
  if (!reticulate::py_numpy_available()) {
    np_ver <- tryCatch(reticulate::import("numpy", convert = TRUE)$`__version__`,
                       error = function(e) "unknown")
    stop("reticulate cannot convert numpy arrays (numpy ", np_ver,
         ", reticulate ", as.character(packageVersion("reticulate")), ").\n",
         "  Upgrade reticulate to 1.38 or newer:  install.packages(\"reticulate\")\n",
         "  or, if you must stay on an older reticulate, install numpy < 2 and ",
         "msprime < 1.4 in the conda environment.")
  }

  if (!quiet) {
    msp <- reticulate::import("msprime", convert = TRUE)
    np  <- reticulate::import("numpy",   convert = TRUE)
    cat(sprintf("conda env : %s\n", conda_env))
    cat(sprintf("python    : %s\n", reticulate::py_config()$version))
    cat(sprintf("msprime   : %s\n", msp$`__version__`))
    cat(sprintf("numpy     : %s\n", np$`__version__`))
  }

  invisible(gs)
}


## ---- record the exact software versions used in a run ----------------------
# Writes a small text file next to the results so a batch can always be traced
# back to the software that produced it. Called automatically by run_batch()
# and run_parallel_batch().
write_session_info <- function(path = "session_info.txt") {
  con <- file(path, open = "wt")
  on.exit(close(con))
  writeLines(c(paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               "", "--- R ---"), con)
  capture.output(print(sessionInfo()), file = con)
  writeLines(c("", "--- Python ---"), con)
  ok <- try({
    msp <- reticulate::import("msprime", convert = TRUE)
    tsk <- reticulate::import("tskit",   convert = TRUE)
    np  <- reticulate::import("numpy",   convert = TRUE)
    writeLines(c(paste("python :", reticulate::py_config()$version),
                 paste("msprime:", msp$`__version__`),
                 paste("tskit  :", tsk$`__version__`),
                 paste("numpy  :", np$`__version__`)), con)
  }, silent = TRUE)
  if (inherits(ok, "try-error")) writeLines("python not initialised", con)
  invisible(path)
}
