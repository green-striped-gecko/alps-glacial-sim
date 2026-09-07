# =============================================================================
# scripts/01_run_simulations.R
# =============================================================================
#
# Runs the full simulation design and writes one .rds per replicate.
#
# Run from the repository root:
#     Rscript scripts/01_run_simulations.R
# or interactively, having set the working directory to the repository root:
#     source("scripts/01_run_simulations.R")
#
# Before the first run, tell reticulate where msprime lives, either by editing
# CONDA_ENV below or by setting the environment variable GLACIALSIM_CONDA_ENV.
#
# WARNING: the published design is 90 replicate runs, each simulating both
# thermal ecologies over 100 time points. Expect a couple of hours per run on
# one core. Start with N_REPS = 1 to check the setup end to end.
# =============================================================================

## ---- settings ---------------------------------------------------------------
CONDA_ENV <- Sys.getenv("GLACIALSIM_CONDA_ENV", unset = "msprime-env")
OUT_DIR   <- "results/batch"
N_REPS    <- 30            # 30 in the paper; use 1 or 2 for a trial
N_CORES   <- 30            # set to 1 to fall back to the serial runner

## ---- load the code ----------------------------------------------------------
source("R/00_setup.R")
source("R/01_simulate.R")
source("R/02_design.R")
source("R/03_run_batch.R")
source("R/04_run_parallel.R")

init_python(CONDA_ENV)

## ---- design -----------------------------------------------------------------
design <- make_design(n_reps = N_REPS)
cat(sprintf("Design: %d rows (%d scenarios x %d replicates)\n",
            nrow(design), length(unique(design$scenario)), N_REPS))

## ---- run --------------------------------------------------------------------
if (N_CORES > 1) {
  run_parallel_batch(design,
                     out_dir   = OUT_DIR,
                     n_cores   = N_CORES,
                     conda_env = CONDA_ENV)
} else {
  run_batch(design, out_dir = OUT_DIR)
}
