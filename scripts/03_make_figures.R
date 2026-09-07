# =============================================================================
# scripts/03_make_figures.R
# =============================================================================
#
# Builds the figures from the simulation output.
#
#   Fig 1       landscape occupancy and per-cell Ne through one glacial cycle
#               needs nothing but the model itself - rebuilt in seconds
#   Fig 2 (H1)  thermal ecology: diversity and structure through time
#               needs the .rds replicate files
#   Fig 5 (H2)  terrain: one peak versus three peaks
#   Fig 6 (H3)  river barrier
#               both need results/pw_dv.rds from scripts/02_extract_summaries.R
#
# The supplementary scenario-comparison figures in R/09_fig_supplementary.R are
# not built here; see that file for how to call them.
#
# Run from the repository root:
#     Rscript scripts/03_make_figures.R
# =============================================================================

SIM_DIRS <- "results/batch"
PW_FILE  <- "results/pw_dv.rds"
FIG_DIR  <- "figures"

source("R/00_setup.R")
source("R/01_simulate.R")      # demog_summary_r()
source("R/03_run_batch.R")     # load_batch()
source("R/06_fig_h1.R")
source("R/07_fig_h2_h3.R")
source("R/08_fig1_occupancy.R")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

## ---- Figure 1: occupancy through a glacial cycle ----------------------------
# Pure demography bookkeeping: this needs the Python engine but no simulation
# results and no batch output, and returns in a few seconds.
init_python(quiet = TRUE)

dat <- build_occupancy_data()
plot_glacial_occupancy(dat, out = file.path(FIG_DIR, "Fig1.png"), dpi = 300)

## ---- Figure 2: H1, thermal ecology ------------------------------------------
results <- load_batch(SIM_DIRS)

# ribbon = "ci95" (the default) is what the published figure used: the 2.5th to
# 97.5th percentile of the replicate means at each time point
plot_h1_trajectories(results,
                     ribbon   = "ci95",
                     png_file = file.path(FIG_DIR, "Fig2_H1_thermal.png"))

## ---- Figures 5 and 6: H2 terrain and H3 river -------------------------------
make_h2_h3_figures(pw_file = PW_FILE, out_dir = FIG_DIR)

cat("\nFigures written to", FIG_DIR, "\n")
