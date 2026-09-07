# =============================================================================
# 02_design.R  --  the simulation design matrix
# =============================================================================
#
# One row per replicate run. Three landscape scenarios x 30 replicates = 90
# rows; each row runs both thermal ecologies, so 180 simulated data sets.
#
# All model parameters are identical across scenarios and replicates. The only
# things that differ are the landscape (`scenario`) and the random seed. That
# is deliberate: the contrasts between scenarios are then attributable to the
# landscape alone.
#
# The parameter values below are the ones used for the published simulations.
# design/published_design_matrix.csv is a saved copy of the result, so the
# exact runs can be reproduced without re-deriving them here.
# =============================================================================


## ---- seeds used for the published runs -------------------------------------
# The same 30 seeds are used for every scenario, so replicate r of "1peak" and
# replicate r of "3peak" share a seed and are directly comparable.
published_seeds <- function() {
  c(1111 + (seq_len(20) - 1L) * 1000L,
    1001 + (seq_len(10) - 1L) * 1000L)
}


## ---- build the design ------------------------------------------------------
#
# n_reps    : replicates per scenario (30 in the paper)
# scenarios : which landscapes to run
# seeds     : one seed per replicate; defaults to the published seeds when
#             n_reps <= 30, otherwise a regular ladder 1001, 2001, 3001, ...
#
# Returns a data.frame ready for run_batch() or run_parallel_batch().
make_design <- function(n_reps    = 30,
                        scenarios = c("1peak", "3peak", "3peak_river"),
                        seeds     = NULL) {

  ## parameters shared by every run -------------------------------------------
  shared <- list(
    n                           = 13,      # 13 x 13 lattice, 169 demes
    cycles                      = 5,       # 5 x 100 ky = 500 ky simulated
    cycle_len                   = 100000,
    warm_len                    = 20000,   # fast warming limb of each cycle
    Ne_good                     = 500,     # deme size at suitability 1
    Ne_bad                      = 50,      # deme size at suitability 0
    m_connected                 = 0.005,
    m_isolated                  = 0.0005,
    seq_len                     = 2e6,     # 2 Mb
    mu                          = 1e-8,
    samples_per_cell_present    = 20,
    samples_per_cell_timeseries = 20,
    ts_step                     = 5000,    # statistics every 5 ky
    do_gl_timeseries            = FALSE,   # genotype time series not retained
    do_pi                       = TRUE
  )

  if (is.null(seeds)) {
    ps <- published_seeds()
    seeds <- if (n_reps <= length(ps)) ps[seq_len(n_reps)]
             else 1001L + (seq_len(n_reps) - 1L) * 1000L
  }
  stopifnot(length(seeds) == n_reps)

  rows <- list()
  for (sc in scenarios) {
    for (r in seq_len(n_reps)) {
      rows[[length(rows) + 1]] <- c(
        list(scenario = sc,
             run      = r,
             seed     = as.integer(seeds[r]),
             file_id  = sprintf("out_%s_run%02d", sc, r)),
        shared)
    }
  }

  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  df$do_gl_timeseries <- as.logical(df$do_gl_timeseries)
  df$do_pi            <- as.logical(df$do_pi)
  rownames(df) <- NULL
  df
}


## ---- turn one design row into arguments for run_both() ---------------------
# Used by both the serial and the parallel runner, so that the two cannot drift
# apart.
design_row_args <- function(row) {
  list(
    scenario                    = as.character(row$scenario),
    n                           = as.integer(row$n),
    cycles                      = as.integer(row$cycles),
    cycle_len                   = as.integer(row$cycle_len),
    warm_len                    = as.integer(row$warm_len),
    Ne_good                     = as.integer(row$Ne_good),
    Ne_bad                      = as.integer(row$Ne_bad),
    m_connected                 = as.numeric(row$m_connected),
    m_isolated                  = as.numeric(row$m_isolated),
    seq_len                     = as.numeric(row$seq_len),
    mu                          = as.numeric(row$mu),
    seed                        = as.integer(row$seed),
    samples_per_cell_present    = as.integer(row$samples_per_cell_present),
    samples_per_cell_timeseries = as.integer(row$samples_per_cell_timeseries),
    ts_step                     = as.integer(row$ts_step),
    do_gl_timeseries            = as.logical(row$do_gl_timeseries),
    do_pi                       = as.logical(row$do_pi)
  )
}
