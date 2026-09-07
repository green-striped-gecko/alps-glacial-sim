# =============================================================================
# 01_simulate.R  --  R wrappers around python/glacial_sim.py
# =============================================================================
#
# These functions convert the output of the Python engine into the R objects
# used downstream: genlight objects (adegenet) for genotype data and plain data
# frames for the time series of summary statistics.
#
# Requires 00_setup.R to have been sourced and init_python() to have been run,
# which creates the module object `gs`.
#
# The scenario ("1peak", "3peak", "3peak_river") is an ordinary argument here.
# Nothing needs to be re-sourced between scenarios.
# =============================================================================

suppressPackageStartupMessages(library(adegenet))

SCENARIOS <- c("1peak", "3peak", "3peak_river")
SPECIES   <- c("warm", "cold")


## ---- genotype matrix -> genlight -------------------------------------------
# G_py comes back from Python as an individuals x sites integer matrix of
# 0/1/2 counts of the derived allele.
.to_genlight <- function(G_py, ind_names, pop_labels, loc_names) {
  G  <- matrix(as.integer(G_py), nrow = length(ind_names))
  gl <- new("genlight", G)
  indNames(gl) <- ind_names
  pop(gl)      <- as.factor(pop_labels)
  locNames(gl) <- loc_names
  gl
}


## ---- 1. present-day genotypes ----------------------------------------------
# One genlight holding samples_per_cell diploid individuals from every deme
# that is occupied at t = 0. Populations are named p<k>, where k is the
# row-major index of the grid cell (k = row * n + column, both 0-based).
gl_present <- function(species, scenario, n, total_time,
                       samples_per_cell, seq_len, mu, seed,
                       Ne_good, Ne_bad, m_connected, m_isolated,
                       cycle_len, warm_len) {

  res <- gs$sim_present(
    n                = as.integer(n),
    species          = species,
    total_time       = as.integer(total_time),
    samples_per_cell = as.integer(samples_per_cell),
    seq_len          = as.numeric(seq_len),
    mu               = mu,
    seed             = as.integer(seed),
    Ne_good          = as.integer(Ne_good),
    Ne_bad           = as.integer(Ne_bad),
    m_connected      = m_connected,
    m_isolated       = m_isolated,
    cycle_len        = as.integer(cycle_len),
    warm_len         = as.integer(warm_len),
    scenario         = scenario
  )

  .to_genlight(
    G_py       = res$G,
    ind_names  = paste0(species, "_t0_", seq_len(res$n_ind)),
    pop_labels = unlist(res$ind_pop),
    loc_names  = paste0("s", as.integer(res$pos))
  )
}


## ---- 2. genotypes at a series of time points -------------------------------
# Returns a named list, one genlight per sampling time. Switched off in the
# published runs (do_gl_timeseries = FALSE) because the statistics needed for
# the figures come from pi_through_time() instead, at a fraction of the cost.
gl_through_time <- function(species, scenario, n, total_time, sample_times,
                            samples_per_cell, seq_len, mu, seed,
                            Ne_good, Ne_bad, m_connected, m_isolated,
                            cycle_len, warm_len) {

  res_list <- gs$sim_timepoints(
    n                = as.integer(n),
    species          = species,
    total_time       = as.integer(total_time),
    sample_times     = as.integer(sample_times),
    samples_per_cell = as.integer(samples_per_cell),
    seq_len          = as.numeric(seq_len),
    mu               = mu,
    seed             = as.integer(seed),
    Ne_good          = as.integer(Ne_good),
    Ne_bad           = as.integer(Ne_bad),
    m_connected      = m_connected,
    m_isolated       = m_isolated,
    cycle_len        = as.integer(cycle_len),
    warm_len         = as.integer(warm_len),
    scenario         = scenario
  )

  blocks <- lapply(res_list, function(b)
    list(t     = as.integer(b$t),
         pop   = as.character(b$pop),
         G     = b$G,
         pos   = as.integer(b$pos),
         n_ind = as.integer(b$n_ind)))

  times_unique <- unique(sapply(blocks, `[[`, "t"))

  gl_list <- lapply(times_unique, function(tt) {
    idx   <- which(sapply(blocks, `[[`, "t") == tt)
    G_all <- do.call(rbind, lapply(blocks[idx], `[[`, "G"))
    pops  <- unlist(lapply(blocks[idx], function(b) rep(b$pop, b$n_ind)))
    pos   <- blocks[[idx[1]]]$pos
    list(time = tt,
         gl   = .to_genlight(G_all,
                             paste0(species, "_t", tt, "_", seq_len(nrow(G_all))),
                             pops, paste0("s", pos)))
  })
  names(gl_list) <- paste0("t", times_unique)
  gl_list
}


## ---- 3. summary statistics through time ------------------------------------
# One row per sampling time. Columns:
#   time      generations before present
#   pi        mean per-deme nucleotide diversity, per base pair
#   tajima_d  mean per-deme Tajima's D
#   theta_w   mean per-deme Watterson's theta, in segregating sites
#   fst_mean  mean Hudson Fst over four-neighbour (adjacent) deme pairs
#   he_mean   mean per-deme expected heterozygosity over segregating sites
pi_through_time <- function(species, scenario, n, total_time, sample_times,
                            samples_per_cell, seq_len, mu, seed,
                            Ne_good, Ne_bad, m_connected, m_isolated,
                            cycle_len, warm_len) {

  res <- gs$sim_pi(
    n                = as.integer(n),
    species          = species,
    total_time       = as.integer(total_time),
    sample_times     = as.integer(sample_times),
    samples_per_cell = as.integer(samples_per_cell),
    seq_len          = as.numeric(seq_len),
    mu               = mu,
    seed             = as.integer(seed),
    Ne_good          = as.integer(Ne_good),
    Ne_bad           = as.integer(Ne_bad),
    m_connected      = m_connected,
    m_isolated       = m_isolated,
    cycle_len        = as.integer(cycle_len),
    warm_len         = as.integer(warm_len),
    scenario         = scenario
  )

  data.frame(time     = as.integer(res$times),
             pi       = as.numeric(res$pi),
             tajima_d = as.numeric(res$tajima_d),
             theta_w  = as.numeric(res$theta_w),
             fst_mean = as.numeric(res$fst_mean),
             he_mean  = as.numeric(res$he_mean))
}


## ---- 4. demography bookkeeping ---------------------------------------------
# No coalescent simulation: just the climate curve, the number of occupied
# demes, the summed and per-deme Ne over time, and the elevation grid.
demog_summary_r <- function(species, scenario, n, total_time, times,
                            Ne_good, Ne_bad, cycle_len, warm_len) {

  res <- gs$demog_summary(
    n          = as.integer(n),
    species    = species,
    total_time = as.integer(total_time),
    times      = as.integer(times),
    Ne_good    = as.integer(Ne_good),
    Ne_bad     = as.integer(Ne_bad),
    cycle_len  = as.integer(cycle_len),
    warm_len   = as.integer(warm_len),
    scenario   = scenario
  )

  list(time    = as.integer(res$times),
       climate = as.numeric(res$climate),
       Ne      = as.numeric(res$Ne_global),
       n_occ   = as.integer(res$n_occ),
       Ne_cell = as.matrix(res$Ne_cell),
       elev    = as.numeric(res$elev))
}


## ---- 5. one complete replicate for one thermal ecology ---------------------
#
# Arguments (defaults are the values used for the published simulations):
#   species        "warm" or "cold"
#   scenario       "1peak", "3peak" or "3peak_river"
#   n              grid side length; 13 gives a 13 x 13 = 169-deme lattice
#   cycles         number of 100 ky glacial cycles simulated
#   cycle_len      length of one climate cycle, generations
#   warm_len       length of the fast warming limb within a cycle
#   Ne_good        deme size at suitability 1
#   Ne_bad         deme size at suitability 0 (refugium size)
#   m_connected    migration rate between occupied neighbouring demes
#   m_isolated     migration rate across a barrier or between isolated refugia
#   seq_len        simulated sequence length, base pairs
#   mu             mutation rate per base pair per generation
#   seed           master seed; the three simulation stages use
#                  seed, seed + 100 and seed + 200
#   samples_per_cell_present     diploid individuals sampled per deme at t = 0
#   samples_per_cell_timeseries  diploid individuals sampled per deme per time point
#   ts_step        spacing of the time series, generations
#   do_gl_timeseries  keep full genotypes at every time point (large; FALSE in
#                     the published runs)
#   do_pi          compute the summary statistic time series (TRUE)
#
# Returns a list with the demography summary, the present-day genlight, the
# optional genotype time series, and the statistics time series.
run_simulation <- function(
    species                     = c("warm", "cold"),
    scenario                    = c("1peak", "3peak", "3peak_river"),
    n                           = 13,
    cycles                      = 5,
    cycle_len                   = 100000,
    warm_len                    = 20000,
    Ne_good                     = 500,
    Ne_bad                      = 50,
    m_connected                 = 0.005,
    m_isolated                  = 0.0005,
    seq_len                     = 2e6,
    mu                          = 1e-8,
    seed                        = 42,
    samples_per_cell_present    = 20,
    samples_per_cell_timeseries = 20,
    ts_step                     = 5000,
    do_gl_timeseries            = FALSE,
    do_pi                       = TRUE) {

  species  <- match.arg(species)
  scenario <- match.arg(scenario)

  total_time <- as.integer(cycles * cycle_len)
  times      <- as.integer(seq(0, total_time - ts_step, by = ts_step))

  cat(sprintf("[%s | %s] n=%d  cycles=%d  Ne_good=%d  Ne_bad=%d  seed=%d\n",
              scenario, species, n, cycles, Ne_good, Ne_bad, seed))

  dem <- demog_summary_r(species, scenario, n, total_time, times,
                         Ne_good, Ne_bad, cycle_len, warm_len)

  cat(sprintf("[%s | %s] present-day genotypes...\n", scenario, species))
  gl_now <- gl_present(species, scenario, n, total_time,
                       samples_per_cell_present, seq_len, mu, seed,
                       Ne_good, Ne_bad, m_connected, m_isolated,
                       cycle_len, warm_len)

  gl_ts <- NULL
  if (do_gl_timeseries) {
    cat(sprintf("[%s | %s] genotypes through time...\n", scenario, species))
    gl_ts <- gl_through_time(species, scenario, n, total_time, times,
                             samples_per_cell_timeseries, seq_len, mu, seed + 100,
                             Ne_good, Ne_bad, m_connected, m_isolated,
                             cycle_len, warm_len)
  }

  pi_df <- NULL
  if (do_pi) {
    cat(sprintf("[%s | %s] diversity statistics through time...\n",
                scenario, species))
    pi_df <- pi_through_time(species, scenario, n, total_time, times,
                             samples_per_cell_timeseries, seq_len, mu, seed + 200,
                             Ne_good, Ne_bad, m_connected, m_isolated,
                             cycle_len, warm_len)
  }

  list(species = species, scenario = scenario, n = n, cycles = cycles,
       total_time = total_time, times = times,
       demog = dem, gl_present = gl_now,
       gl_timeseries = gl_ts, pi = pi_df)
}


## ---- 6. both thermal ecologies on the same landscape -----------------------
# The warm- and cold-adapted taxa share the landscape, the climate history and
# the seed, so a replicate is a matched pair.
run_both <- function(scenario = c("1peak", "3peak", "3peak_river"), ...) {
  scenario <- match.arg(scenario)
  cat("=== WARM-ADAPTED ===\n"); flush.console()
  warm <- run_simulation(species = "warm", scenario = scenario, ...)
  cat("=== COLD-ADAPTED ===\n"); flush.console()
  cold <- run_simulation(species = "cold", scenario = scenario, ...)
  list(warm = warm, cold = cold)
}
