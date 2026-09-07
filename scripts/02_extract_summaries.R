# =============================================================================
# scripts/02_extract_summaries.R
# =============================================================================
#
# Recomputes present-day all-pairs Fst and per-deme diversity from the .rds
# files written by scripts/01_run_simulations.R, and saves them as one compact
# object (pw_dv.rds) that the H2/H3 figure code reads.
#
# This step needs no Python: it works on the genlight objects already stored in
# the results. It is memory-hungry rather than slow, because each replicate
# holds a full present-day genotype matrix; the extraction deliberately keeps
# only one ecology in memory at a time.
#
# Run from the repository root:
#     Rscript scripts/02_extract_summaries.R
# =============================================================================

SIM_DIRS <- "results/batch"          # one or more directories of out_*.rds
OUT_FILE <- "results/pw_dv.rds"

source("R/05_extract_fst.R")

if (!dir.exists(dirname(OUT_FILE))) dir.create(dirname(OUT_FILE), recursive = TRUE)

extract_tables(sim_dirs = SIM_DIRS, out_file = OUT_FILE)
