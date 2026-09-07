# =============================================================================
# scripts/00_check_install.R
# =============================================================================
#
# Quick check that everything needed is present and talking to everything else.
# Run this first on a new machine:
#     Rscript scripts/00_check_install.R
#
# It reports the versions of every dependency, then runs one very small
# simulation (5 x 5 grid, one cycle, three time points) end to end. That takes
# well under a minute. If it finishes, the full pipeline will run.
# =============================================================================

CONDA_ENV <- Sys.getenv("GLACIALSIM_CONDA_ENV", unset = "msprime-env")

cat("== R packages ==\n")
need <- c("reticulate", "adegenet", "dplyr", "tibble", "tidyr",
          "ggplot2", "patchwork", "mgcv", "scales", "parallel")
for (p in need) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-12s %s\n", p,
              if (ok) as.character(packageVersion(p)) else "MISSING"))
}
if (!all(sapply(need, requireNamespace, quietly = TRUE)))
  stop("Install the missing R packages listed above.")

cat("\n== Python ==\n")
source("R/00_setup.R")
init_python(CONDA_ENV)

cat("\n== Test simulation ==\n")
source("R/01_simulate.R")

t0 <- proc.time()["elapsed"]
test <- run_simulation(species = "cold", scenario = "3peak_river",
                       n = 5, cycles = 1, ts_step = 25000,
                       Ne_good = 100, Ne_bad = 20,
                       seq_len = 2e5, seed = 1,
                       samples_per_cell_present = 4,
                       samples_per_cell_timeseries = 4)

cat(sprintf("\nfinished in %.1f s\n", proc.time()["elapsed"] - t0))
cat("present-day genlight:", nInd(test$gl_present), "individuals,",
    nLoc(test$gl_present), "loci,",
    nlevels(pop(test$gl_present)), "demes\n")
print(test$pi)

cat("\nInstallation looks good.\n")
