# =============================================================================
# 06_fig_h1.R  --  H1 figure: diversity and structure through time
# =============================================================================
#
# Two stacked panels sharing one time axis, with a climate strip on top:
#   A  mean per-deme nucleotide diversity (pi)
#   B  mean Fst between four-neighbour deme pairs
#
# Every curve is a mean over replicates; the ribbon is the spread across
# replicates. Because the ribbon is taken over independent replicate runs, and
# each replicate contributes one value per time point, the interval reflects
# run-to-run variance rather than variance among demes within a run.
#
# Encoding
#   colour family : thermal ecology  (warm = orange, cold = blue)
#   tint          : landscape scenario (light -> dark: 1peak, 3peak, +river)
#   line type     : landscape scenario (solid, dashed, dotted)
#
# The warm-adapted curve for "3peak_river" is dropped by default. During the
# glacial phases in which the river exists the warm-adapted taxon is confined to
# the lowland corner refugia, so the barrier never touches it and the curve is
# identical to "3peak" by construction rather than by result.
#
# Base graphics only; no ggplot dependency.
# =============================================================================

COL_WARM <- "#D2691E"
COL_COLD <- "#4682B4"

# three tints per ecology, one per scenario
.SP_TINT <- list(
  warm = c("1peak" = "#F0A868", "3peak" = "#D2691E", "3peak_river" = "#8A3C08"),
  cold = c("1peak" = "#8FBBDA", "3peak" = "#4682B4", "3peak_river" = "#1F4E79"))
.SC_LTY   <- c("1peak" = 1, "3peak" = 2, "3peak_river" = 3)
.SC_LABEL <- c("1peak" = "1-peak", "3peak" = "3-peak",
               "3peak_river" = "3-peak + river")


## ---- statistic matrix: time points x replicates -----------------------------
.mat <- function(results, scenario, species, stat) {
  keep <- grep(paste0("^out_", scenario, "_run"), names(results))
  sapply(results[keep], function(x) x[[species]]$pi[[stat]])
}

## ---- replicate mean and interval at each time point -------------------------
# ribbon = "ci95"   2.5th to 97.5th percentile across replicates
#          "sd"     mean +/- one standard deviation across replicates
#          "minmax" full range across replicates
.stats <- function(results, sc, sp, stat, ribbon = "ci95") {
  m <- .mat(results, sc, sp, stat)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  mu <- rowMeans(m, na.rm = TRUE)
  if (ncol(m) < 2) return(list(mean = mu, lo = NULL, hi = NULL, n = ncol(m)))
  if (ribbon == "minmax") {
    lo <- apply(m, 1, min, na.rm = TRUE); hi <- apply(m, 1, max, na.rm = TRUE)
  } else if (ribbon == "sd") {
    s <- apply(m, 1, sd, na.rm = TRUE); lo <- mu - s; hi <- mu + s
  } else {
    lo <- apply(m, 1, quantile, 0.025, na.rm = TRUE)
    hi <- apply(m, 1, quantile, 0.975, na.rm = TRUE)
  }
  list(mean = mu, lo = lo, hi = hi, n = ncol(m))
}


## ---- the figure -------------------------------------------------------------
#
# results   : named list from load_batch()
# scenarios : which landscapes to draw, in legend order
# ribbon    : "ci95", "sd" or "minmax" (see .stats above)
# pi_max    : upper limit of panel A; curves above it are simply clipped
# fst_max   : upper limit of panel B, NULL to take it from the data
# t_max     : oldest time shown, NULL for the whole simulation
# drop_identical_warm : omit warm + river (see header)
# clip_marks : mark with a triangle where a mean curve leaves panel B, so a
#              capped axis is not mistaken for the data ending
# png_file  : write a 300 dpi PNG instead of drawing to the active device
#
# Returns, invisibly, the list of per-curve means and intervals.
plot_h1_trajectories <- function(
    results,
    scenarios = c("1peak", "3peak", "3peak_river"),
    ribbon    = c("ci95", "sd", "minmax"),
    pi_max    = 0.008,
    fst_max   = NULL,
    t_max     = NULL,
    drop_identical_warm = TRUE,
    cycle_len = 100000, warm_len = 20000,
    cex_lab = 1.15, cex_leg = 1.15, cex_ax = 0.9,
    ribbon_alpha = 0.30, lwd_line = 1.8,
    clip_marks = TRUE,
    png_file = NULL, width = 8, height = 8.5) {

  ribbon <- match.arg(ribbon)
  if (!is.null(png_file))
    png(png_file, width = width, height = height, units = "in", res = 300)
  op <- par(no.readonly = TRUE)
  on.exit({par(op); if (!is.null(png_file)) dev.off()})

  S <- list()
  for (sc in scenarios) for (sp in c("warm", "cold"))
    S[[paste(sc, sp)]] <- list(
      pi  = .stats(results, sc, sp, "pi",       ribbon),
      fst = .stats(results, sc, sp, "fst_mean", ribbon))

  nrep <- sapply(scenarios, function(sc)
    length(grep(paste0("^out_", sc, "_run"), names(results))))
  message("replicates per scenario: ",
          paste(sprintf("%s=%d", scenarios, nrep), collapse = ", "),
          "; ribbon = ", ribbon)
  if (any(nrep < 2))
    warning("scenarios with <2 replicates draw no ribbon: ",
            paste(scenarios[nrep < 2], collapse = ", "), call. = FALSE)

  first <- results[[grep(paste0("^out_", scenarios[1], "_run"),
                         names(results))[1]]]
  tw <- first$warm$pi$time
  if (is.null(t_max)) t_max <- max(tw)
  xlim <- c(t_max, 0)                       # time runs right to left

  curves <- expand.grid(sc = scenarios, sp = c("warm", "cold"),
                        stringsAsFactors = FALSE)
  if (drop_identical_warm && all(c("3peak", "3peak_river") %in% scenarios))
    curves <- curves[!(curves$sp == "warm" & curves$sc == "3peak_river"), ]
  .col <- function(i) .SP_TINT[[curves$sp[i]]][curves$sc[i]]

  t_fine <- seq(0, t_max, by = 1000)
  clim <- sapply(t_fine, function(t) {
    x <- t %% cycle_len
    if (x <= warm_len) 1 - x / warm_len else
      (x - warm_len) / (cycle_len - warm_len) })

  # glacial phases shaded using the panel's own ylim, not par("usr")
  .glacials <- function(yl) {
    r <- rle(clim < 0.5); p <- cumsum(c(1, r$lengths))
    for (i in seq_along(r$lengths)) if (r$values[i])
      rect(t_fine[p[i]], yl[1] - diff(yl),
           t_fine[min(p[i] + r$lengths[i] - 1, length(t_fine))],
           yl[2] + diff(yl),
           col = adjustcolor("#4682B4", 0.11), border = NA)
  }
  .vlines <- function() abline(v = t_fine[clim == 1 | clim == 0],
                               lty = 3, col = "grey75", lwd = 0.7)

  .draw <- function(key, yl, mark = FALSE) {
    plot(NULL, xlim = xlim, ylim = yl, xaxt = "n", xlab = "", ylab = "",
         las = 1, cex.axis = cex_ax)
    .glacials(yl); .vlines()
    for (i in seq_len(nrow(curves))) {
      s <- S[[paste(curves$sc[i], curves$sp[i])]][[key]]
      if (!is.null(s$lo))
        polygon(c(tw, rev(tw)), c(s$lo, rev(s$hi)),
                col = adjustcolor(.col(i), ribbon_alpha),
                border = adjustcolor(.col(i), 0.45), lwd = 0.4)
    }
    for (i in seq_len(nrow(curves))) {
      s <- S[[paste(curves$sc[i], curves$sp[i])]][[key]]
      lines(tw, s$mean, lwd = lwd_line, col = .col(i),
            lty = .SC_LTY[curves$sc[i]])
    }
    if (mark) for (i in seq_len(nrow(curves))) {
      s   <- S[[paste(curves$sc[i], curves$sp[i])]][[key]]
      out <- which(s$mean > yl[2])
      if (length(out)) {
        grp <- cumsum(c(1, diff(out) != 1))       # one mark per excursion
        xs  <- tapply(tw[out], grp, function(v) mean(range(v)))
        points(xs, rep(yl[2], length(xs)), pch = 17, cex = 0.7,
               col = .col(i), xpd = NA)
      }
    }
    box()
  }

  layout(matrix(1:3), heights = c(0.7, 2.6, 2.6))
  par(mar = c(0, 6.0, 2.2, 2.5), oma = c(5.2, 0, 0, 0))

  ## climate strip
  plot(t_fine, clim, type = "l", lwd = 2, col = "grey25", xlim = xlim,
       ylim = c(-.05, 1.2), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  .glacials(c(-.05, 1.2)); lines(t_fine, clim, lwd = 2, col = "grey25"); .vlines()
  mtext("Climate", side = 2, line = 1.6, cex = cex_lab * 0.8)
  u <- par("usr")
  text(u[1], 1.10, " Interglacial", col = "grey40", adj = c(0, .5),
       cex = cex_leg * 0.8, font = 3)
  text(u[1], -0.02, " Glacial", col = "grey40", adj = c(0, 0),
       cex = cex_leg * 0.8, font = 3)

  ## panel A: nucleotide diversity
  par(mar = c(0, 6.0, 1.4, 2.5))
  .draw("pi", c(0, pi_max))
  mtext(expression("Nucleotide diversity " * pi), side = 2, line = 4.2,
        cex = cex_lab)
  mtext("A", side = 3, line = 0.1, adj = -0.10, font = 2, cex = cex_lab * 1.1)

  # bespoke key: one swatch per scenario tint, one row per ecology
  .tint_key <- function(y_npc, species, label) {
    seg <- 0.038; gap <- 0.008; x0 <- 0.030
    y <- grconvertY(y_npc, "npc", "user")
    for (k in seq_along(scenarios)) {
      xa <- grconvertX(x0 + (k - 1) * (seg + gap), "npc", "user")
      xb <- grconvertX(x0 + (k - 1) * (seg + gap) + seg, "npc", "user")
      segments(xa, y, xb, y, lwd = 3.6, lend = 1,
               col = .SP_TINT[[species]][scenarios[k]])
    }
    text(grconvertX(x0 + length(scenarios) * (seg + gap) + 0.012, "npc", "user"),
         y, label, adj = c(0, 0.5), cex = cex_leg, xpd = NA)
  }
  .tint_key(0.94, "warm", "warm-adapted")
  .tint_key(0.86, "cold", "cold-adapted")
  legend("topright", legend = .SC_LABEL[scenarios], col = "grey35",
         lty = .SC_LTY[scenarios], lwd = 2.4, bty = "n",
         cex = cex_leg, x.intersp = 0.8, seg.len = 2.4)

  ## panel B: Fst between adjacent demes
  par(mar = c(0, 6.0, 1.4, 2.5))
  fv <- unlist(lapply(S, function(s) c(s$fst$mean, s$fst$hi)))
  fv <- fv[is.finite(fv)]
  if (is.null(fst_max)) fst_max <- min(1, max(fv) * 1.06)
  .draw("fst", c(0, fst_max), mark = clip_marks)
  mtext(expression("Mean pairwise " * F[ST]), side = 2, line = 4.2, cex = cex_lab)
  mtext("B", side = 3, line = 0.1, adj = -0.10, font = 2, cex = cex_lab * 1.1)

  ## shared x axis
  span  <- t_max / 1000
  sp_ky <- if (span <= 120) 20 else if (span <= 300) 50 else 100
  tk    <- seq(0, floor(span / sp_ky) * sp_ky, by = sp_ky) * 1000
  axis(1, at = tk, labels = ifelse(tk == 0, "present", format(tk / 1000)),
       cex.axis = cex_ax)
  mtext("Thousands of years before present", side = 1, line = 3.0, cex = cex_lab)

  invisible(S)
}
