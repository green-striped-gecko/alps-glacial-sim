# =============================================================================
# 07_fig_h2_h3.R  --  H2 (terrain) and H3 (river barrier) figures
# =============================================================================
#
# Input is the list(pw, dv) produced by R/05_extract_fst.R, i.e. present-day
# all-pairs Fst and per-deme diversity for every replicate.
#
# Figure 5 (H2, one peak versus three peaks)
#   A  distribution of pairwise Fst by scenario and ecology
#   B  distribution of per-deme pi by scenario and ecology
#   C  Fst against geographic distance, one panel per scenario
#
# Figure 6 (H3, three peaks with and without the river)
#   A  Fst against distance, split by whether the pair spans the river
#   B  Fst excess of river-spanning pairs over the isolation-by-distance
#      expectation fitted on non-spanning pairs
#   C  per-deme pi by scenario and ecology
#
# Two modelling choices are worth stating explicitly.
#
# Isolation by distance is concave: Fst saturates as distance grows. A straight
# line fitted mostly on short-distance pairs is too steep, overpredicts long
# pairs, and can run above Fst = 1. A shrinkage spline is used instead, with
# the basis dimension set from the number of distinct distances so that sparse
# strata degrade gracefully to something close to a straight line.
#
# The barrier effect is measured as a residual, not as a raw difference. The
# isolation-by-distance curve is fitted on non-spanning pairs only and then
# used to predict the spanning pairs; the residual is the excess structure
# attributable to the river once distance is accounted for. Spanning pairs
# further apart than any non-spanning pair are dropped, because there the
# spline would be extrapolating.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(mgcv)
})


## ---- isolation-by-distance smoother ----------------------------------------
ibd_k <- function(x, cap = 6) max(3, min(cap, length(unique(x)) - 1))

ibd_form <- function(x, cap = 6)
  as.formula(sprintf("y ~ s(x, bs = 'cs', k = %d)", ibd_k(x, cap)))

ibd_fit <- function(d, cap = 6) {
  k <- ibd_k(d$dist, cap)
  if (nrow(d) < 10 || k < 3) return(lm(fst ~ dist, data = d))
  mgcv::gam(fst ~ s(dist, bs = "cs", k = k), data = d)
}


## ---- shared styling ---------------------------------------------------------
EC <- c(`Cold-adapted` = "#2C6E9B", `Warm-adapted` = "#C1553B")

h2h3_theme <- function()
  theme_bw(9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position  = "bottom",
        plot.tag         = element_text(face = "bold", size = 10))


## ---- main entry point -------------------------------------------------------
#
# pw_file  : path to the list(pw, dv) written by extract_tables()
# out_dir  : where the PDF and PNG versions of the figures are written
# prefix   : file name stem, e.g. "Fig5_H2_terrain"
#
# Returns, invisibly, the summary tables printed at the end.
make_h2_h3_figures <- function(pw_file = "pw_dv.rds",
                               out_dir = ".",
                               fig5_name = "Fig5_H2_terrain",
                               fig6_name = "Fig6_H3_river") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  theme_set(h2h3_theme())

  d  <- readRDS(pw_file)
  pw <- d$pw
  dv <- d$dv

  ## ---- river position ------------------------------------------------------
  # The simulation puts the river on Python grid row (n %/% 2) - 1, 0-based,
  # which is row n %/% 2 in the 1-based coordinates of these tables. Demes
  # sitting *on* the river cannot be assigned to a bank, so their pairs get
  # crosses = NA rather than being forced to one side.
  GRID_N <- unique(pw$grid_n)
  stopifnot(length(GRID_N) == 1)
  RIVER_ROW <- GRID_N %/% 2
  side <- function(rr) sign(rr - RIVER_ROW)          # -1 / 0 / +1

  pw <- pw %>%
    mutate(crosses = ifelse(side(row_i) == 0 | side(row_j) == 0, NA,
                            side(row_i) != side(row_j)),
           ecology = factor(ecology, c("cold", "warm"),
                            c("Cold-adapted", "Warm-adapted")))

  ## =========================================================================
  ## FIGURE 5 : H2, terrain
  ## =========================================================================
  h2 <- pw %>%
    filter(scenario %in% c("1peak", "3peak")) %>%
    mutate(scenario = factor(scenario, c("1peak", "3peak"),
                             c("One peak", "Three peaks")))

  h2d <- dv %>%
    filter(scenario %in% c("1peak", "3peak")) %>%
    mutate(scenario = factor(scenario, c("1peak", "3peak"),
                             c("One peak", "Three peaks")),
           ecology  = factor(ecology, c("cold", "warm"),
                             c("Cold-adapted", "Warm-adapted")))

  pA <- ggplot(h2, aes(scenario, fst, fill = ecology)) +
    geom_violin(scale = "width", alpha = .65, linewidth = .25,
                colour = "grey30", position = position_dodge(width = .9)) +
    stat_summary(aes(group = ecology), fun = mean, geom = "point", size = 1.7,
                 shape = 21, colour = "grey15", fill = "white", stroke = .5,
                 position = position_dodge(width = .9)) +
    scale_fill_manual(values = EC, name = NULL) +
    labs(x = NULL, y = expression("pairwise " * F[ST]), tag = "A")

  pB <- ggplot(h2d, aes(scenario, pi, fill = ecology)) +
    geom_violin(scale = "width", alpha = .65, linewidth = .25,
                colour = "grey30", position = position_dodge(width = .9)) +
    stat_summary(aes(group = ecology), fun = mean, geom = "point", size = 1.7,
                 shape = 21, colour = "grey15", fill = "white", stroke = .5,
                 position = position_dodge(width = .9)) +
    scale_fill_manual(values = EC, name = NULL) +
    scale_y_continuous(labels = function(v) v * 1000) +
    labs(x = NULL,
         y = expression("per-deme " * pi * " (per site, " %*% 10^-3 * ")"),
         tag = "B")

  pC <- ggplot(h2, aes(dist, fst, colour = ecology)) +
    geom_point(alpha = .06, size = .5) +
    geom_smooth(method = "gam", formula = ibd_form(h2$dist), se = TRUE,
                linewidth = .6) +
    facet_wrap(~scenario) +
    scale_colour_manual(values = EC, name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "geographic distance (cells)", y = expression(F[ST]), tag = "C")

  fig5 <- (pA | pB) / pC +
    plot_layout(heights = c(1, 1.1), guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(file.path(out_dir, paste0(fig5_name, ".pdf")), fig5,
         width = 174, height = 150, units = "mm")
  ggsave(file.path(out_dir, paste0(fig5_name, ".png")), fig5,
         width = 174, height = 150, units = "mm", dpi = 200)

  ## =========================================================================
  ## FIGURE 6 : H3, river barrier
  ## =========================================================================
  n_on_river <- sum(is.na(pw$crosses[pw$scenario %in% c("3peak", "3peak_river")]))
  message("dropped ", n_on_river, " pairs with a deme on the river row")

  h3 <- pw %>%
    filter(scenario %in% c("3peak", "3peak_river"), !is.na(crosses)) %>%
    mutate(scenario = factor(scenario, c("3peak", "3peak_river"),
                             c("Three peaks", "Three peaks + river")),
           crosses  = factor(crosses, c(FALSE, TRUE), c("no", "yes")))

  pD <- ggplot(h3, aes(dist, fst, colour = crosses)) +
    geom_point(alpha = .10, size = .5) +
    geom_smooth(method = "gam", formula = ibd_form(h3$dist), se = TRUE,
                linewidth = .6) +
    facet_grid(ecology ~ scenario) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(values = c(no = "grey55", yes = "#B03A2E"),
                        name = "pair spans river") +
    labs(x = "geographic distance (cells)", y = expression(F[ST]), tag = "A")

  # residualise against the isolation-by-distance expectation of non-spanning pairs
  h3 <- h3 %>%
    group_by(scenario, ecology) %>%
    group_modify(function(dd, key) {
      b <- filter(dd, crosses == "no")
      if (nrow(b) > 5) {
        dd$resid <- as.numeric(dd$fst - predict(ibd_fit(b), newdata = dd))
        dd$resid[dd$dist > max(b$dist)] <- NA_real_   # no spline support here
      } else dd$resid <- NA_real_
      dd
    }) %>% ungroup()

  pE <- ggplot(h3 %>% filter(crosses == "yes"),
               aes(scenario, resid, fill = ecology)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_violin(scale = "width", alpha = .65, linewidth = .25, colour = "grey30") +
    stat_summary(fun = mean, geom = "point", size = 1.7, shape = 21,
                 colour = "grey15", fill = "white", stroke = .5) +
    facet_wrap(~ecology) +
    scale_fill_manual(values = EC, guide = "none") +
    labs(x = NULL,
         y = expression(F[ST] * " excess in river-spanning pairs"), tag = "B") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  dvr <- dv %>%
    filter(scenario %in% c("3peak", "3peak_river")) %>%
    mutate(ecology  = factor(ecology, c("cold", "warm"),
                             c("Cold-adapted", "Warm-adapted")),
           scenario = factor(scenario, c("3peak", "3peak_river"),
                             c("Three peaks", "Three peaks + river")))

  pF <- ggplot(dvr, aes(scenario, pi, fill = ecology)) +
    geom_violin(scale = "width", alpha = .65, linewidth = .25, colour = "grey30") +
    stat_summary(fun = mean, geom = "point", size = 1.7, shape = 21,
                 colour = "grey15", fill = "white", stroke = .5) +
    facet_wrap(~ecology) +
    scale_fill_manual(values = EC, guide = "none") +
    labs(x = NULL,
         y = expression("per-deme " * pi * " (per site, " %*% 10^-3 * ")"),
         tag = "C") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  fig6 <- pD / (pE | pF) + plot_layout(heights = c(1.3, 1))

  ggsave(file.path(out_dir, paste0(fig6_name, ".pdf")), fig6,
         width = 174, height = 170, units = "mm")
  ggsave(file.path(out_dir, paste0(fig6_name, ".png")), fig6,
         width = 174, height = 170, units = "mm", dpi = 200)

  ## ---- numbers quoted in the text ------------------------------------------
  s_h2 <- h2 %>% group_by(scenario, ecology) %>%
    summarise(mean_fst = round(mean(fst, na.rm = TRUE), 3), .groups = "drop")
  s_h3 <- h3 %>% filter(crosses == "yes") %>%
    group_by(scenario, ecology) %>%
    summarise(excess = round(mean(resid, na.rm = TRUE), 4), n = n(),
              .groups = "drop")

  cat("\n--- H2: mean pairwise Fst ---\n"); print(s_h2)
  cat("\n--- H3: Fst excess of river-spanning pairs over IBD expectation ---\n")
  print(s_h3)

  invisible(list(h2 = s_h2, h3 = s_h3))
}
