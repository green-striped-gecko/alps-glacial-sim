# =============================================================================
# 09_fig_supplementary.R  --  earlier scenario-comparison figures
# =============================================================================
#
# Four exploratory figures that preceded the published Figure 2. They plot all
# six curves (3 landscapes x 2 thermal ecologies) on one axis, with colour for
# ecology and line type for landscape. Figure 2 replaced them with the stacked
# two-panel design in R/06_fig_h1.R, which separates diversity from structure.
#
#   plot_pi_scenarios()          nucleotide diversity through time
#   plot_he_scenarios()          expected heterozygosity through time
#   plot_fst_scenarios()         adjacent-pair Fst through time
#   plot_fst_present_scenarios() present-day pairwise Fst, boxplots
#
# Usage:
#   results <- load_batch("results/batch")
#   plot_pi_scenarios(results, ribbon = "minmax", png_file = "SFig_pi.png")
#
# NOTE ON DEPENDENCIES
# plot_fst_present_scenarios() is the only function in the whole repository
# that needs dartRverse: it calls dartR.base::gl.report.fstat(). The other
# three, and everything else in the pipeline, run on adegenet alone. The
# present-day Fst used for Figures 5 and 6 is computed independently in
# R/05_extract_fst.R with Hudson's estimator, so dartRverse is optional.
#
# Internal helpers and constants are prefixed .sup_ so that sourcing this file
# alongside R/06_fig_h1.R cannot overwrite the styling used by Figure 2.
# =============================================================================

.sup_SP_COL <- c(warm = COL_WARM, cold = COL_COLD)

# ── scenario linetypes: solid=1peak, dashed=3peak, dotted=3peak_river
.sup_SC_LTY <- c("1peak" = 1, "3peak" = 2, "3peak_river" = 3)
.sup_SC_LABEL <- c(
  "1peak"       = "1-peak",
  "3peak"       = "3-peak",
  "3peak_river" = "3-peak river"
)


# ── internal: mean + ribbon bounds for one scenario/species ───
.sup_sc_stats <- function(results, scenario, species, stat = "pi",
                      ribbon = "sd") {
  mat <- extract_pi_matrix(results, scenario, species, stat)
  m   <- rowMeans(mat, na.rm = TRUE)
  if (ribbon == "minmax") {
    lo <- apply(mat, 1, min, na.rm = TRUE)
    hi <- apply(mat, 1, max, na.rm = TRUE)
  } else {
    s  <- apply(mat, 1, sd, na.rm = TRUE)
    lo <- m - s
    hi <- m + s
  }
  list(mean = m, lo = lo, hi = hi, n = ncol(mat))
}


# ── kybp x-axis (right = present, left = deep past) ──────────
.sup_kybp_axis <- function(tw) {
  span_ky <- max(tw) / 1000
  spacing <- if (span_ky <= 150) 20 else if (span_ky <= 400) 50 else 100
  tick_ky <- seq(0, ceiling(span_ky / spacing) * spacing, by = spacing)
  tick_yr <- tick_ky * 1000
  keep    <- tick_yr <= max(tw)          # filter both together
  tick_yr <- tick_yr[keep]
  tick_ky <- tick_ky[keep]
  labs    <- ifelse(tick_ky == 0, "present", format(paste0(tick_ky,",000")))
  axis(1, at = tick_yr, labels = labs, cex.axis = 0.82)
  mtext(" Years Before Present", side = 1, line = 3.0, cex = 0.88)
}


plot_pi_scenarios <- function(
    results,
    scenarios = c("1peak", "3peak", "3peak_river"),
    ribbon    = c("sd", "minmax"),
    log_pi    = FALSE,
    png_file  = NULL,
    width     = 9,
    height    = 6.5
) {
  ribbon <- match.arg(ribbon)
  
  if (!is.null(png_file))
    png(png_file, width = width, height = height, units = "in", res = 300)
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2), heights = c(1, 3.5))
  par(mar = c(0, 5.5, 2, 2), oma = c(5, 0, 0, 0))
  
  # ── pull stats ────────────────────────────────────────────────
  stats <- list()
  for (sc in scenarios) {
    stats[[sc]] <- list(
      warm = .sup_sc_stats(results, sc, "warm", ribbon = ribbon),
      cold = .sup_sc_stats(results, sc, "cold", ribbon = ribbon)
    )
  }
  
  # ── shared time / climate from first replicate ────────────────
  first <- results[[grep(paste0("out_", scenarios[1], "_run"),
                         names(results))[1]]]
  tw    <- first$warm$pi$time
  tc    <- first$warm$demog$time
  clim  <- first$warm$demog$climate
  xlim  <- rev(range(tw))
  
  # ── shared y range ────────────────────────────────────────────
  all_vals <- unlist(lapply(stats, function(s)
    c(s$warm$hi, s$cold$hi, s$warm$mean, s$cold$mean)))
  
  if (log_pi) {
    pos_vals <- all_vals[all_vals > 0]
    ylim     <- c(max(min(pos_vals, na.rm = TRUE) * 0.5, 1e-10),
                  max(all_vals, na.rm = TRUE) * 2)
    log_arg  <- "y"
  } else {
    ylim    <- c(0, max(all_vals, na.rm = TRUE) * 1.18)
    log_arg <- ""
  }
  
  # ── helpers ───────────────────────────────────────────────────
  # Use a fine time grid (1k steps) for shading so the clim=0.5
  # crossing is accurate — coarse demography steps miss the exact
  # boundary (e.g. clim=0.5 is at t=10k, not t=20k).
  t_fine   <- seq(0, max(tc), by = 1000)
  cycle_len <- 100000L; warm_len <- 20000L
  clim_fine <- sapply(t_fine, function(t) {
    x <- t %% cycle_len
    if (x <= warm_len) 1.0 - x / warm_len else (x - warm_len) / (cycle_len - warm_len)
  })
  
  .glacials <- function() {
    usr  <- par("usr")
    in_g <- clim_fine < 0.5
    runs <- rle(in_g); pos <- cumsum(c(1, runs$lengths))
    for (i in seq_along(runs$lengths)) {
      if (!runs$values[i]) next
      t0 <- t_fine[pos[i]]
      t1 <- t_fine[min(pos[i] + runs$lengths[i] - 1, length(t_fine))]
      rect(t0, usr[3], t1, usr[4],
           col = adjustcolor("#4682B4", 0.11), border = NA)
    }
  }
  .vlines <- function() {
    # dotted lines at every clim=1 (warm peak) and clim=0 (LGM)
    abline(v = t_fine[clim_fine == 1], lty = 3, col = "grey70", lwd = 0.8)
    abline(v = t_fine[clim_fine == 0], lty = 3, col = "grey70", lwd = 0.8)
  }
  .ribbon_poly <- function(x, lo, hi, col, alpha = 0.18) {
    lo[is.na(lo)] <- 0; hi[is.na(hi)] <- 0
    polygon(c(x, rev(x)), c(lo, rev(hi)),
            col = adjustcolor(col, alpha), border = NA)
  }
  
  # ── Panel 1: climate ─────────────────────────────────────────
  plot(t_fine, clim_fine,
       type = "l", lwd = 2, col = "grey25",
       xlim = xlim, ylim = c(-0.05, 1.20),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Climate")
  .glacials(); lines(t_fine, clim_fine, lwd = 2, col = "grey25")
  abline(h = c(0, 1), lty = 2, col = "grey70", lwd = 0.8)
  .vlines()
  usr <- par("usr")
  text(usr[1], 1.10, " Interglacial", col = "grey40", adj = c(0, .5), cex = 0.68, font = 3)
  text(usr[1], -0.03, " Glacial",      col = "grey40", adj = c(0, -.5), cex = 0.68, font = 3)
  
  ribbon_lbl <- if (ribbon == "sd") "ribbon = \u00b11 SD" else "ribbon = min\u2013max"
  mtext(bquote(bold("Nucleotide Diversity  ") * pi * bold("(t)") ~
                 "  " * .(ribbon_lbl)),
        side = 3, line = 0.5, cex = 0.92, font = 2, adj = 0.5)
  
  # ── Panel 2: all 6 π curves ───────────────────────────────────
  par(mar = c(0, 5.5, 0, 2))
  
  plot(NULL, xlim = xlim, ylim = ylim, log = log_arg,
       xaxt = "n", xlab = "", ylab = expression("Nucleotide Diversity " * pi))
  .glacials(); .vlines()
  
  # ribbons first (behind lines)
  for (sc in scenarios) {
    .ribbon_poly(tw, stats[[sc]]$cold$lo, stats[[sc]]$cold$hi, .sup_SP_COL["cold"])
    .ribbon_poly(tw, stats[[sc]]$warm$lo, stats[[sc]]$warm$hi, .sup_SP_COL["warm"])
  }
  
  # mean lines: cold beneath warm per scenario
  for (sc in scenarios) {
    lty <- .sup_SC_LTY[sc]
    lines(tw, stats[[sc]]$cold$mean, lwd = 2.0, col = .sup_SP_COL["cold"], lty = lty)
    lines(tw, stats[[sc]]$warm$mean, lwd = 2.0, col = .sup_SP_COL["warm"], lty = lty)
  }
  
  # ── kybp x axis ──────────────────────────────────────────────
  .sup_kybp_axis(tw)
  
  # ── legend ───────────────────────────────────────────────────
  sc_ok <- scenarios[scenarios %in% names(.sup_SC_LTY)]
  legend(
    "topleft",
    legend = c("warm-adapted", "cold-adapted", .sup_SC_LABEL[sc_ok]),
    col    = c(.sup_SP_COL["warm"], .sup_SP_COL["cold"], rep("grey30", length(sc_ok))),
    lty    = c(1, 1, .sup_SC_LTY[sc_ok]),
    lwd    = c(2.5, 2.5, rep(2.0, length(sc_ok))),
    bty = "n", cex = 0.76, x.intersp = 0.7
  )
  
  if (!is.null(png_file)) dev.off()
  invisible(stats)
}




plot_he_scenarios <- function(
    results,
    scenarios = c("1peak", "3peak", "3peak_river"),
    ribbon    = c("sd", "minmax"),
    log_he    = FALSE,
    png_file  = NULL,
    width     = 9,
    height    = 6.5
) {
  ribbon <- match.arg(ribbon)
  
  if (!is.null(png_file))
    png(png_file, width = width, height = height, units = "in", res = 150)
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2), heights = c(1, 3.5))
  par(mar = c(0, 5.5, 2, 2), oma = c(5, 0, 0, 0))
  
  # ── pull He stats ─────────────────────────────────────────────
  stats <- list()
  for (sc in scenarios) {
    stats[[sc]] <- list(
      warm = .sup_sc_stats(results, sc, "warm", stat = "he_mean", ribbon = ribbon),
      cold = .sup_sc_stats(results, sc, "cold", stat = "he_mean", ribbon = ribbon)
    )
  }
  
  # ── shared time / climate from first replicate ────────────────
  first <- results[[grep(paste0("out_", scenarios[1], "_run"),
                         names(results))[1]]]
  tw    <- first$warm$pi$time
  tc    <- first$warm$demog$time
  clim  <- first$warm$demog$climate
  xlim  <- rev(range(tw))
  
  # ── shared y range ────────────────────────────────────────────
  all_vals <- unlist(lapply(stats, function(s)
    c(s$warm$hi, s$cold$hi, s$warm$mean, s$cold$mean)))
  all_vals <- all_vals[is.finite(all_vals)]
  
  if (log_he) {
    pos_vals <- all_vals[all_vals > 0]
    ylim     <- c(max(min(pos_vals, na.rm = TRUE) * 0.5, 1e-10),
                  max(all_vals, na.rm = TRUE) * 2)
    log_arg  <- "y"
  } else {
    ylim    <- c(0.2, max(all_vals, na.rm = TRUE) * 1.18)
    log_arg <- ""
  }
  
  # ── fine climate grid for accurate shading ────────────────────
  t_fine    <- seq(0, max(tc), by = 1000)
  cycle_len <- 100000L; warm_len <- 20000L
  clim_fine <- sapply(t_fine, function(t) {
    x <- t %% cycle_len
    if (x <= warm_len) 1.0 - x / warm_len else (x - warm_len) / (cycle_len - warm_len)
  })
  
  .glacials <- function() {
    usr  <- par("usr")
    in_g <- clim_fine < 0.5
    runs <- rle(in_g); pos <- cumsum(c(1, runs$lengths))
    for (i in seq_along(runs$lengths)) {
      if (!runs$values[i]) next
      t0 <- t_fine[pos[i]]
      t1 <- t_fine[min(pos[i] + runs$lengths[i] - 1, length(t_fine))]
      rect(t0, usr[3], t1, usr[4],
           col = adjustcolor("#4682B4", 0.11), border = NA)
    }
  }
  .vlines <- function() {
    abline(v = t_fine[clim_fine == 1], lty = 3, col = "grey70", lwd = 0.8)
    abline(v = t_fine[clim_fine == 0], lty = 3, col = "grey70", lwd = 0.8)
  }
  .ribbon_poly <- function(x, lo, hi, col, alpha = 0.18) {
    lo[is.na(lo)] <- 0; hi[is.na(hi)] <- 0
    polygon(c(x, rev(x)), c(lo, rev(hi)),
            col = adjustcolor(col, alpha), border = NA)
  }
  
  # ── Panel 1: climate ─────────────────────────────────────────
  plot(t_fine, clim_fine,
       type = "l", lwd = 2, col = "grey25",
       xlim = xlim, ylim = c(-0.05, 1.20),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Climate")
  .glacials(); lines(t_fine, clim_fine, lwd = 2, col = "grey25")
  abline(h = c(0, 1), lty = 2, col = "grey70", lwd = 0.8)
  .vlines()
  usr <- par("usr")
  text(usr[1], 1.10, " Interglacial", col = "grey40", adj = c(0, .5), cex = 0.68, font = 3)
  text(usr[1], -0.03, " Glacial",      col = "grey40", adj = c(0, -0.5), cex = 0.68, font = 3)
  
  ribbon_lbl <- if (ribbon == "sd") "ribbon = \u00b11 SD" else "ribbon = min\u2013max"
  mtext(bquote(bold("Expected heterozygosity  ") * H[e] * bold("(t)") ~
                 "  " * .(ribbon_lbl)),
        side = 3, line = 0.5, cex = 0.92, font = 2, adj = 0.5)
  
  # ── Panel 2: all 6 He curves ──────────────────────────────────
  par(mar = c(0, 5.5, 0, 2))
  
  plot(NULL, xlim = xlim, ylim = ylim, log = log_arg,
       xaxt = "n", xlab = "", ylab = expression("Expected Heterozygosity " * H[e]))
  .glacials(); .vlines()
  
  # ribbons first (behind lines)
  for (sc in scenarios) {
    .ribbon_poly(tw, stats[[sc]]$cold$lo, stats[[sc]]$cold$hi, .sup_SP_COL["cold"])
    .ribbon_poly(tw, stats[[sc]]$warm$lo, stats[[sc]]$warm$hi, .sup_SP_COL["warm"])
  }
  
  # mean lines: cold beneath warm per scenario
  for (sc in scenarios) {
    lty <- .sup_SC_LTY[sc]
    lines(tw, stats[[sc]]$cold$mean, lwd = 2.0, col = .sup_SP_COL["cold"], lty = lty)
    lines(tw, stats[[sc]]$warm$mean, lwd = 2.0, col = .sup_SP_COL["warm"], lty = lty)
  }
  
  # ── kybp x axis ──────────────────────────────────────────────
  .sup_kybp_axis(tw)
  
  # ── legend ───────────────────────────────────────────────────
  sc_ok <- scenarios[scenarios %in% names(.sup_SC_LTY)]
  legend(
    "topleft",
    legend = c("warm-adapted", "cold-adapted", .sup_SC_LABEL[sc_ok]),
    col    = c(.sup_SP_COL["warm"], .sup_SP_COL["cold"], rep("grey30", length(sc_ok))),
    lty    = c(1, 1, .sup_SC_LTY[sc_ok]),
    lwd    = c(2.5, 2.5, rep(2.0, length(sc_ok))),
    bty = "n", cex = 0.76, x.intersp = 0.7
  )
  
  if (!is.null(png_file)) dev.off()
  invisible(stats)
}


# ============================================================
# plot_fst_scenarios.R
# Present-day pairwise Fst boxplots: 3 scenarios × 2 species
# Pools pairwise Fst values across ALL replicates per scenario.
#
# Requires: dartRverse (gl.fst.pop), plot_simple.R (COL_WARM, COL_COLD),
#           run_batch.R (load_batch)
#
# Usage:
#   results <- load_batch("batch_output")
#   plot_fst_scenarios(results)
#   plot_fst_scenarios(results, png_file = "fst_scenarios.png")
#
# Layout: 6 boxplots, paired warm+cold per scenario.
#   Colour: warm = COL_WARM, cold = COL_COLD.
#   Scenario label centred above each warm+cold pair.
#   Dashed vertical line separating scenarios.
#   Each box pools all pairwise Fst values across all replicates.
# ============================================================


plot_fst_present_scenarios <- function(
    results,
    scenarios = c("1peak", "3peak", "3peak_river"),
    png_file  = NULL,
    width     = 8,
    height    = 6
) {
  
  if (!requireNamespace("dartRverse", quietly = TRUE))
    stop("dartRverse required: install.packages('dartRverse' )")
  
  if (!is.null(png_file))
    png(png_file, width = width, height = height, units = "in", res = 300)
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  
  # ── extract pairwise Fst from one genlight object ────────────
  .fst_vals <- function(gl) {
    fst_mat <- suppressMessages(
      dartR.base::gl.report.fstat(gl[,], verbose = 0)$Stat_matrices$Fst
    )
    vals <- fst_mat[upper.tri(fst_mat)]
    vals <- vals[!is.na(vals)]
    pmax(vals, 0)
  }
  
  # ── pool Fst across all replicates per scenario × species ─────
  .pool_fst <- function(scenario, species) {
    keys <- grep(paste0("out_", scenario, "_run"), names(results), value = TRUE)
    if (length(keys) == 0)
      stop("No replicates found for scenario: ", scenario)
    cat(sprintf("  %s / %s: %d replicates\n", scenario, species, length(keys)))
    vals <- unlist(lapply(keys, function(k) {
      tryCatch(.fst_vals(results[[k]][[species]]$gl_present),
               error = function(e) { message("  skipping ", k, ": ", e$message); NULL })
    }))
    vals[!is.na(vals)]
  }
  
  # ── compute all 6 pooled distributions ───────────────────────
  cat("Computing Fst across replicates...\n")
  sc_labels_nice <- c("1peak"       = "1-peak",
                      "3peak"       = "3-peak",
                      "3peak_river" = "3-peak river")
  
  fst  <- list()
  keys <- character()
  for (sc in scenarios) {
    for (sp in c("warm", "cold")) {
      key       <- paste0(sc, "\n", sp)
      fst[[key]] <- .pool_fst(sc, sp)
      keys      <- c(keys, key)
    }
  }
  
  n_reps <- length(grep(paste0("out_", scenarios[1], "_run"), names(results)))
  
  # ── colours ───────────────────────────────────────────────────
  box_cols  <- rep(c(adjustcolor(COL_WARM, 0.35),
                     adjustcolor(COL_COLD,  0.35)), length(scenarios))
  bord_cols <- rep(c(COL_WARM, COL_COLD), length(scenarios))
  
  # ── y range ──────────────────────────────────────────────────
  ymax <- max(unlist(fst), na.rm = TRUE) * 1.25
  ymax <- max(ymax, 0.20)
  
  # ── plot ─────────────────────────────────────────────────────
  par(mar = c(4, 5, 4.5, 2))
  
  bp <- boxplot(fst,
                col      = box_cols,
                border   = bord_cols,
                names    = rep("", length(fst)),
                ylab     = expression(F[ST] ~ "(pairwise, all occupied cell pairs)"),
                ylim     = c(0, ymax),
                outline  = FALSE,
                whisklty = 1,
                las      = 1,
                main     = "")
  
  # ── jittered points ───────────────────────────────────────────
  set.seed(1)
  for (i in seq_along(fst)) {
    jx <- jitter(rep(i, length(fst[[i]])), amount = 0.12)
    points(jx, fst[[i]],
           pch = 19, cex = 0.30,
           col = adjustcolor(bord_cols[i], 0.25))
  }
  
  # ── reference lines ───────────────────────────────────────────
  abline(h = c(0.05, 0.15, 0.25), lty = 3, col = "grey75")
  text(par("usr")[2], c(0.05, 0.15, 0.25),
       c("low", "mod", "high"),
       adj = c(1, -0.35), cex = 0.60, col = "grey55", font = 3, xpd = TRUE)
  
  # ── x-axis species labels ─────────────────────────────────────
  axis(1, at = seq_along(fst),
       labels = rep(c("warm", "cold"), length(scenarios)),
       tick = FALSE, line = 0, cex.axis = 0.82)
  
  # ── scenario labels centred above each pair ───────────────────
  sc_at <- seq(1.5, by = 2, length.out = length(scenarios))
  for (i in seq_along(scenarios)) {
    mtext(sc_labels_nice[scenarios[i]], side = 3, at = sc_at[i],
          line = 0.4, cex = 0.88, font = 2)
  }
  
  # ── dashed separators between scenarios ──────────────────────
  sep_at <- seq(2.5, by = 2, length.out = length(scenarios) - 1)
  abline(v = sep_at, lty = 2, col = "grey50", lwd = 1.2)
  
  # ── median annotations ────────────────────────────────────────
  for (i in seq_along(fst)) {
    med <- median(fst[[i]], na.rm = TRUE)
    text(i, med, sprintf("%.3f", med),
         adj = c(0.5, -0.5), cex = 0.62, font = 2,
         col = bord_cols[i])
  }
  
  # ── legend ────────────────────────────────────────────────────
  legend("topright",
         legend = c("warm-adapted", "cold-adapted"),
         fill   = c(adjustcolor(COL_WARM, 0.35), adjustcolor(COL_COLD, 0.35)),
         border = c(COL_WARM, COL_COLD),
         bty = "n", cex = 0.80)
  
  mtext(expression(bold("Present-day pairwise  ") * F[ST] ~ "(t = 0)"),
        side = 3, line = 2.8, cex = 1.0, font = 2)
  mtext(sprintf("(pooled across %d replicates per scenario)", n_reps),
        side = 3, line = 1.6, cex = 0.78, col = "grey40")
  
  if (!is.null(png_file)) dev.off()
  
  # ── summary table ─────────────────────────────────────────────
  summary_df <- data.frame(
    scenario   = rep(scenarios, each = 2),
    species    = rep(c("warm", "cold"), length(scenarios)),
    n_reps     = n_reps,
    n_pairs    = sapply(fst, length),
    fst_median = round(sapply(fst, median, na.rm = TRUE), 4),
    fst_mean   = round(sapply(fst, mean,   na.rm = TRUE), 4),
    fst_max    = round(sapply(fst, max,    na.rm = TRUE), 4),
    row.names  = NULL
  )
  cat("\n"); print(summary_df)
  invisible(summary_df)
}

### new function

plot_fst_scenarios <- function(
    results,
    scenarios = c("1peak", "3peak", "3peak_river"),
    ribbon    = c("sd", "minmax"),
    log_fst   = FALSE,
    png_file  = NULL,
    width     = 9,
    height    = 6.5
) {
  ribbon <- match.arg(ribbon)
  
  if (!is.null(png_file))
    png(png_file, width = width, height = height, units = "in", res = 300)
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2), heights = c(1, 3.5))
  par(mar = c(0, 5.5, 2, 2), oma = c(5, 0, 0, 0))
  
  # ── pull Fst stats ────────────────────────────────────────────
  stats <- list()
  for (sc in scenarios) {
    stats[[sc]] <- list(
      warm = .sup_sc_stats(results, sc, "warm", stat = "fst_mean", ribbon = ribbon),
      cold = .sup_sc_stats(results, sc, "cold", stat = "fst_mean", ribbon = ribbon)
    )
  }
  
  # ── shared time / climate from first replicate ────────────────
  first <- results[[grep(paste0("out_", scenarios[1], "_run"),
                         names(results))[1]]]
  tw    <- first$warm$pi$time
  tc    <- first$warm$demog$time
  xlim  <- rev(range(tw))
  
  # ── shared y range ────────────────────────────────────────────
  all_vals <- unlist(lapply(stats, function(s)
    c(s$warm$hi, s$cold$hi, s$warm$mean, s$cold$mean)))
  all_vals <- all_vals[is.finite(all_vals)]
  
  if (log_fst) {
    pos_vals <- all_vals[all_vals > 0]
    ylim     <- c(max(min(pos_vals, na.rm = TRUE) * 0.5, 1e-6),
                  max(all_vals, na.rm = TRUE) * 2)
    log_arg  <- "y"
  } else {
    ylim    <- c(0, max(all_vals, na.rm = TRUE) * 1.18)
    log_arg <- ""
  }
  
  # ── fine climate grid ─────────────────────────────────────────
  t_fine    <- seq(0, max(tc), by = 1000)
  cycle_len <- 100000L; warm_len <- 20000L
  clim_fine <- sapply(t_fine, function(t) {
    x <- t %% cycle_len
    if (x <= warm_len) 1.0 - x / warm_len else (x - warm_len) / (cycle_len - warm_len)
  })
  
  .glacials <- function() {
    usr  <- par("usr")
    in_g <- clim_fine < 0.5
    runs <- rle(in_g); pos <- cumsum(c(1, runs$lengths))
    for (i in seq_along(runs$lengths)) {
      if (!runs$values[i]) next
      t0 <- t_fine[pos[i]]
      t1 <- t_fine[min(pos[i] + runs$lengths[i] - 1, length(t_fine))]
      rect(t0, usr[3], t1, usr[4],
           col = adjustcolor("#4682B4", 0.11), border = NA)
    }
  }
  .vlines <- function() {
    abline(v = t_fine[clim_fine == 1], lty = 3, col = "grey70", lwd = 0.8)
    abline(v = t_fine[clim_fine == 0], lty = 3, col = "grey70", lwd = 0.8)
  }
  .ribbon_poly <- function(x, lo, hi, col, alpha = 0.18) {
    lo[is.na(lo)] <- 0; hi[is.na(hi)] <- 0
    polygon(c(x, rev(x)), c(lo, rev(hi)),
            col = adjustcolor(col, alpha), border = NA)
  }
  
  # ── Panel 1: climate ──────────────────────────────────────────
  plot(t_fine, clim_fine,
       type = "l", lwd = 2, col = "grey25",
       xlim = xlim, ylim = c(-0.05, 1.20),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "Climate")
  .glacials(); lines(t_fine, clim_fine, lwd = 2, col = "grey25")
  abline(h = c(0, 1), lty = 2, col = "grey70", lwd = 0.8)
  .vlines()
  usr <- par("usr")
  text(usr[1], 1.10, " Interglacial", col = "grey40", adj = c(0, .5), cex = 0.68, font = 3)
  text(usr[1], -0.03, " Glacial",      col = "grey40", adj = c(0, -0.5), cex = 0.68, font = 3)
  
  ribbon_lbl <- if (ribbon == "sd") "ribbon = \u00b11 SD" else "ribbon = min\u2013max"
  mtext(bquote(bold("Mean pairwise  ") * F[ST] * bold("(t)") ~
                 "  " * .(ribbon_lbl)),
        side = 3, line = 0.5, cex = 0.92, font = 2, adj = 0.5)
  
  # ── Panel 2: all 6 Fst curves ─────────────────────────────────
  par(mar = c(0, 5.5, 0, 2))
  
  plot(NULL, xlim = xlim, ylim = ylim, log = log_arg,
       xaxt = "n", xlab = "", ylab = expression(F[ST]))
  .glacials(); .vlines()
  
  # reference lines
  abline(h = c(0.05, 0.15, 0.25), lty = 3, col = "grey80", lwd = 0.8)
  
  # ribbons first
  for (sc in scenarios) {
    .ribbon_poly(tw, stats[[sc]]$cold$lo, stats[[sc]]$cold$hi, .sup_SP_COL["cold"])
    .ribbon_poly(tw, stats[[sc]]$warm$lo, stats[[sc]]$warm$hi, .sup_SP_COL["warm"])
  }
  
  # mean lines
  for (sc in scenarios) {
    lty <- .sup_SC_LTY[sc]
    lines(tw, stats[[sc]]$cold$mean, lwd = 2.0, col = .sup_SP_COL["cold"], lty = lty)
    lines(tw, stats[[sc]]$warm$mean, lwd = 2.0, col = .sup_SP_COL["warm"], lty = lty)
  }
  
  # ── kybp x axis ──────────────────────────────────────────────
  .sup_kybp_axis(tw)
  
  # ── legend ───────────────────────────────────────────────────
  sc_ok <- scenarios[scenarios %in% names(.sup_SC_LTY)]
  legend(
    "topleft",
    legend = c("warm-adapted", "cold-adapted", .sup_SC_LABEL[sc_ok]),
    col    = c(.sup_SP_COL["warm"], .sup_SP_COL["cold"], rep("grey30", length(sc_ok))),
    lty    = c(1, 1, .sup_SC_LTY[sc_ok]),
    lwd    = c(2.5, 2.5, rep(2.0, length(sc_ok))),
    bty = "n", cex = 0.76, x.intersp = 0.7
  )
  
  if (!is.null(png_file)) dev.off()
  invisible(stats)
}