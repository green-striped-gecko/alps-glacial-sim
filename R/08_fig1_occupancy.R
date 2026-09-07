# =============================================================================
# 08_fig1_occupancy.R  --  Figure 1: landscape occupancy through a glacial cycle
# =============================================================================
#
# A climate strip across the top, then four rows of per-cell Ne heatmaps at
# 10 ky intervals through one 100 ky cycle, running 100k (left) to the present
# (right):
#
#   1-peak          warm-adapted, then cold-adapted
#   3-peak(+river)  warm-adapted, then cold-adapted, with the river drawn
#
# Warm and cold use separate red and blue colour ramps sharing one Ne range;
# unoccupied cells (Ne <= ne_min) are left white, so the figure shows both how
# large the demes are and where the taxon is at all.
#
# The river is drawn as a horizontal line between the northern massif and the
# two southern ones: dotted and thin while it is inactive, solid and thick
# while the climate is glacial (c < 0.5) and the barrier exists.
#
# Everything here is demography bookkeeping. No coalescent simulation is
# involved, so the figure can be rebuilt in seconds with
# build_occupancy_data() below; the batch results are not needed.
#
# Requires: ggplot2, patchwork, dplyr, scales
# =============================================================================


## ---- build the input from the model -----------------------------------------
#
# Returns the four-element list that plot_glacial_occupancy() expects, in the
# order it expects: 1-peak warm, 1-peak cold, 3-peak+river warm, 3-peak+river
# cold. Each element carries $n, $times and $demog (climate, Ne_cell, elev).
#
# Defaults match the published simulations. Needs init_python() to have run,
# because demog_summary_r() calls the Python engine, but it is pure
# bookkeeping and returns in a couple of seconds.
build_occupancy_data <- function(n         = 13,
                                 cycles    = 5,
                                 cycle_len = 100000,
                                 warm_len  = 20000,
                                 Ne_good   = 500,
                                 Ne_bad    = 50,
                                 ts_step   = 5000) {

  total_time <- as.integer(cycles * cycle_len)
  times      <- as.integer(seq(0, total_time - ts_step, by = ts_step))

  combos <- list(list(sp = "warm", sc = "1peak"),
                 list(sp = "cold", sc = "1peak"),
                 list(sp = "warm", sc = "3peak_river"),
                 list(sp = "cold", sc = "3peak_river"))

  lapply(combos, function(k) {
    dem <- demog_summary_r(k$sp, k$sc, n, total_time, times,
                           Ne_good, Ne_bad, cycle_len, warm_len)
    list(species = k$sp, scenario = k$sc, n = n,
         cycles = cycles, total_time = total_time, times = times,
         demog = dem)
  })
}


## ---- the figure --------------------------------------------------------------
#
# dat        four-element list from build_occupancy_data(), or the equivalent
#            slice of four run_simulation() outputs
# out        output path; the device is chosen from the extension
#            (.pdf, .png or .svg)
# width,
# height     inches
# dpi        resolution for raster output
# ne_min     cells at or below this Ne count as unoccupied and are drawn white
# snap_every generations between snapshot columns
# rotate     "auto" orients the 3-peak landscape so the lone peak is north;
#            or give 0, 90, 180 or 270 degrees explicitly
# base_size  global font size
plot_glacial_occupancy <- function(
    dat,
    out        = "glacial_occupancy.pdf",
    width      = 14,
    height     = 9,
    dpi        = 300,
    ne_min     = 2,
    snap_every = 10000,
    rotate     = "auto",
    base_size  = 11
) {

  # ---- 0. packages ----------------------------------------------------------
  for (pkg in c("ggplot2", "patchwork", "dplyr", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf(
        "Package '%s' is required. Install with: install.packages('%s')",
        pkg, pkg))
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  WARM_COL <- "#C94B2A"
  COLD_COL <- "#2A6EA6"

  # ---- 1. output device from the file extension -----------------------------
  ext <- tolower(tools::file_ext(out))
  if (!ext %in% c("pdf", "png", "svg"))
    stop("Output format must be one of: .pdf, .png, .svg  (got '.", ext, "')")

  device <- switch(ext,
                   pdf = cairo_pdf,
                   png = function(filename, ...)
                     png(filename, ..., units = "in", res = dpi),
                   svg = svg)

  # ---- 2. snapshot indices --------------------------------------------------
  # Simulation times are 0, 5000, ..., 495000 (five 100 ky cycles). One cycle
  # runs t = 0 (c = 1, interglacial) -> t = 20000 (c = 0, glacial maximum) ->
  # t = 100000 (c = 1 again). Simulation time IS years before present, so
  # panels run 100k on the left to the present on the right. t = 100000 has the
  # same climate state as t = 0, so it reuses the t = 0 index.
  all_times <- dat[[1]]$times
  snap_t    <- seq(0, 100000, by = snap_every)
  snap_idx  <- match(snap_t %% 100000, all_times)
  kybp      <- snap_t / 1000
  ord       <- order(-kybp)
  snap_t <- snap_t[ord]; snap_idx <- snap_idx[ord]; kybp <- kybp[ord]
  clim_snap <- dat[[1]]$demog$climate[snap_idx]
  N         <- dat[[1]]$n

  # ---- 2b. landscape geometry and river position ----------------------------
  # CAREFUL: make_df fills df$row from the matrix COLUMN index and df$col from
  # the matrix ROW index - they are transposed relative to their names. With
  # aes(x = row, y = col) the plotted vertical axis is therefore the matrix
  # ROW. So the peaks must split 1-against-2 across matrix rows, with the lone
  # peak at the smaller row, and the river is drawn at that row on the y axis.
  .rot <- function(m, deg) switch(
    as.character(deg),
    "0" = m, "90" = t(m[N:1, ]), "180" = m[N:1, N:1], "270" = t(m)[N:1, ],
    stop("rotate must be one of 0, 90, 180, 270"))

  .peak_mat <- function(elem) {
    z <- elem$demog$elev
    if (is.null(z)) z <- apply(elem$demog$Ne_cell, 2, max, na.rm = TRUE)
    matrix(z, nrow = N, ncol = N, byrow = TRUE)      # as in make_df
  }

  # how the three summits split across matrix rows (the plotted vertical axis)
  .split <- function(m) {
    pk <- which(m >= max(m) * 0.999, arr.ind = TRUE)
    if (nrow(pk) < 3) return(NULL)
    tb <- table(pk[, "row"])
    if (length(tb) != 2 || !any(tb == 1)) return(NULL)
    list(lone = as.numeric(names(tb)[tb == 1]),
         pair = as.numeric(names(tb)[tb != 1]))
  }

  m3 <- NULL
  for (elem in dat) {
    mm <- .peak_mat(elem)
    if (!is.null(.split(mm)) || !is.null(.split(t(mm)))) { m3 <- mm; break }
  }

  river_at <- NULL
  if (is.null(m3)) {
    if (identical(rotate, "auto")) rotate <- 0
    message("geometry: no 3-peak landscape found; no river drawn")
  } else {
    if (identical(rotate, "auto")) {
      rotate <- NA
      for (d in c(0, 90, 180, 270)) {
        sp <- .split(.rot(m3, d))
        if (!is.null(sp) && sp$lone < sp$pair) { rotate <- d; break }
      }
      if (is.na(rotate)) {
        rotate <- 0
        warning("could not orient the landscape automatically; using 0 deg",
                call. = FALSE)
      }
    }
    sp <- .split(.rot(m3, rotate))
    if (is.null(sp)) {
      message(sprintf(paste("geometry: rotation %d deg leaves no north-south",
                            "peak split; no river drawn"), rotate))
    } else {
      river_at <- mean(c(sp$lone, sp$pair))
      message(sprintf(paste("geometry: rotate %d deg -> lone peak row %g",
                            "(north), paired peaks row %g (south),",
                            "river at y = %.1f"),
                      rotate, sp$lone, sp$pair, river_at))
    }
  }

  # ---- 3. tidy data frame for one row of panels -----------------------------
  make_df <- function(elem, deg = 0) {
    rows <- lapply(seq_along(snap_idx), function(s) {
      ne_vec <- elem$demog$Ne_cell[snap_idx[s], ]
      ne_mat <- matrix(ne_vec, nrow = N, ncol = N, byrow = TRUE)
      ne_mat <- .rot(ne_mat, deg)
      ne_mat[ne_mat <= ne_min] <- NA          # unoccupied -> white
      data.frame(
        row  = rep(seq_len(N), N),
        col  = rep(seq_len(N), each = N),
        Ne   = as.vector(t(ne_mat)),          # row-major order
        kybp = kybp[s],
        clim = clim_snap[s])
    })
    df      <- do.call(rbind, rows)
    df$kybp <- factor(df$kybp, levels = kybp)  # left (100k) -> right (0k)
    df
  }

  # only the 3-peak panels need rotating; the 1-peak landscape is symmetric
  dfs <- lapply(dat, make_df, deg = rotate)

  # ---- 4. colour scales -----------------------------------------------------
  ne_max <- max(sapply(dat, function(e) max(e$demog$Ne_cell, na.rm = TRUE)))

  make_fill <- function(is_warm) {
    scale_fill_gradientn(
      colours = if (is_warm) c("#fee0d2", "#fc9272", "#de2d26", "#a50f15")
                else         c("#deebf7", "#9ecae1", "#3182bd", "#08306b"),
      limits   = c(ne_min + 0.01, ne_max),
      na.value = "white",
      guide    = "none")
  }

  # ---- 5. shared heatmap theme ----------------------------------------------
  theme_hm <- function(show_strip = FALSE) {
    theme_minimal(base_size = base_size) +
      theme(
        panel.spacing = unit(0.4, "mm"),
        strip.text    = if (show_strip)
          element_text(size = base_size * 0.92, colour = "#333333",
                       lineheight = 1.1, margin = margin(b = 2))
          else element_blank(),
        strip.clip   = "off",
        axis.text    = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(colour = "#cccccc", fill = NA,
                                    linewidth = 0.3),
        plot.margin  = margin(1, 4, 1, 8, "mm"))
  }

  # ---- 5b. river geometry ---------------------------------------------------
  # Panels are built with byrow = TRUE, so the barrier lies along a grid ROW,
  # separating the single peak from the pair. It is only active while the
  # climate is glacial (c < 0.5); dotted otherwise.
  df_river <- if (is.null(river_at)) NULL else data.frame(
    kybp   = factor(kybp, levels = kybp),
    riv    = river_at,
    active = factor(clim_snap < 0.5, levels = c("FALSE", "TRUE")))

  # ---- 6. heatmap factory ---------------------------------------------------
  make_heatmap <- function(df, is_warm, row_label, show_strip = FALSE,
                           show_river = FALSE) {
    strip_labs <- setNames(sprintf("%dk", kybp), as.character(kybp))
    ggplot(df, aes(x = row, y = col, fill = Ne)) +
      geom_tile(colour = NA) +
      make_fill(is_warm) +
      facet_wrap(~ kybp, nrow = 1, labeller = labeller(kybp = strip_labs)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_reverse(expand = c(0, 0)) +
      { if (show_river && !is.null(df_river)) list(
          geom_segment(data = df_river, inherit.aes = FALSE,   # white underlay
                       aes(x = 0.5, xend = N + 0.5, y = riv, yend = riv),
                       colour = "white", linewidth = 1.6),
          geom_segment(data = df_river, inherit.aes = FALSE,
                       aes(x = 0.5, xend = N + 0.5, y = riv, yend = riv,
                           linewidth = active, linetype = active),
                       colour = "#0B7285", lineend = "butt"),
          scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.2),
                                 guide = "none"),
          scale_linetype_manual(values = c(`FALSE` = "22", `TRUE` = "solid"),
                                guide = "none")) else NULL } +
      labs(y = row_label) +
      theme_hm(show_strip) +
      theme(axis.title.y = element_text(
        size = base_size * 1.45, face = "bold", angle = 90, vjust = 0.5,
        lineheight = 0.95,
        colour = if (is_warm) WARM_COL else COLD_COL,
        margin = margin(r = 3)))
  }

  # ---- 7. climate panel -----------------------------------------------------
  t_fine  <- seq(0, 100000, by = 200)
  clim_fn <- function(t) {
    x <- t %% 100000
    ifelse(x <= 20000, 1 - x / 20000, (x - 20000) / 80000)
  }
  df_c <- data.frame(kybp = t_fine / 1000, clim = clim_fn(t_fine))

  p_clim <- ggplot(df_c, aes(x = kybp, y = clim)) +
    geom_ribbon(aes(ymin = 0, ymax = pmin(clim, 0.5)),
                fill = COLD_COL, alpha = 0.35) +
    geom_ribbon(aes(ymin = pmin(clim, 0.5), ymax = clim),
                fill = WARM_COL, alpha = 0.35) +
    geom_line(colour = "#222222", linewidth = 0.75) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "#888888", linewidth = 0.4) +
    geom_vline(xintercept = kybp, colour = "#cccccc",
               linewidth = 0.3, alpha = 0.9) +
    # glacial maximum: c = 0 at t = 20000, i.e. 20 kybp on this axis
    geom_vline(xintercept = 20, colour = COLD_COL,
               linewidth = 0.75, linetype = "dotted") +
    annotate("text", x = 21, y = 1.10, label = "LGM",
             colour = COLD_COL, size = base_size * 0.62, fontface = "bold",
             hjust = 0) +
    geom_vline(xintercept = 0, colour = WARM_COL,
               linewidth = 0.75, linetype = "dotted") +
    annotate("text", x = 1.2, y = 1.10, label = "Present\n(interglacial)",
             colour = WARM_COL, size = base_size * 0.60, fontface = "bold",
             hjust = 1, lineheight = 0.9) +
    scale_x_reverse(breaks = kybp, labels = paste0(kybp, "k"),
                    limits = c(max(kybp) + snap_every / 2000,
                               min(kybp) - snap_every / 2000),
                    expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1),
                       labels = c("0\n(cold)", "0.5", "1\n(warm)"),
                       limits = c(-0.04, 1.22)) +
    labs(x = "kyr before present", y = "Climate\nindex c") +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x      = element_text(size = base_size * 1.55,
                                      colour = "#333333"),
      axis.text.y      = element_text(size = base_size * 1.20),
      axis.title.x     = element_text(size = base_size * 1.30),
      axis.title.y     = element_text(size = base_size * 1.45, angle = 90),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#f0f0f0", linewidth = 0.3),
      plot.margin      = margin(4, 4, 1, 8, "mm"))

  # ---- 8. heatmap rows ------------------------------------------------------
  p_w1  <- make_heatmap(dfs[[1]], TRUE,  "Warm")
  p_c1  <- make_heatmap(dfs[[2]], FALSE, "Cold")
  p_w3r <- make_heatmap(dfs[[3]], TRUE,  "Warm", show_river = TRUE)
  p_c3r <- make_heatmap(dfs[[4]], FALSE, "Cold", show_river = TRUE)

  # ---- 9. assemble ----------------------------------------------------------
  .side <- function(txt) wrap_elements(full = grid::textGrob(
    txt, rot = 90, gp = grid::gpar(fontsize = base_size * 1.5,
                                   fontface = "bold", col = "#333333")))
  lab1 <- .side("1-peak"); lab2 <- .side("3-peak(+river)")

  design <- "
AB
CD
CE
FG
FH
"
  final <- plot_spacer() + p_clim + lab1 + p_w1 + p_c1 + lab2 + p_w3r + p_c3r +
    plot_layout(design = design, widths = c(1, 26),
                heights = c(1.9, 1, 1, 1, 1)) +
    plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2, "mm")))

  # ---- 10. save -------------------------------------------------------------
  ggsave(filename = out, plot = final, device = device,
         width = width, height = height, bg = "white")
  message(sprintf("Saved [%s] -> %s", toupper(ext), out))
  invisible(final)
}
