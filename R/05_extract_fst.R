# =============================================================================
# 05_extract_fst.R  --  present-day Fst and per-deme diversity
# =============================================================================
#
# The time series written by the simulation (the $pi data frame) contains an
# Fst averaged over *four-neighbour deme pairs only*. The tests of H2 (terrain)
# and H3 (river barrier) need *all* pairs of demes at the present day, together
# with the geographic distance between them, so that isolation by distance can
# be separated from the barrier effect.
#
# This file recomputes Fst from the present-day genlight objects. Two tables
# come out:
#
#   pw : one row per pair of demes per replicate per ecology
#        run, scenario, ecology, deme_i, deme_j, row_i, col_i, row_j, col_j,
#        dist (Euclidean, in grid cells), fst, grid_n
#
#   dv : one row per deme per replicate per ecology
#        run, scenario, ecology, deme, row, col, pi (per base pair),
#        n_snp, seq_len, grid_n
#
# Estimator
# ---------
# Hudson's Fst as a ratio of averages, which is the low-bias choice when demes
# differ in sample size and many SNPs are near-monomorphic (Bhatia et al. 2013).
# For demes i and j with allele frequencies p and sample sizes N:
#
#   num_ij = sum_l (p_il - p_jl)^2 - c_i - c_j ,  c_i = sum_l p_il(1-p_il)/(N_il-1)
#   den_ij = sum_l [ p_il(1-p_jl) + p_jl(1-p_il) ]
#   Fst_ij = num_ij / den_ij
#
# Both numerator and denominator expand into P %*% t(P) plus per-deme row sums,
# so the whole matrix of pairwise values costs one matrix multiplication.
#
# Grid coordinates
# ----------------
# Populations are named p<k> with k = row * n + column, both 0-based, matching
# cell_idx() in python/glacial_sim.py. Here they are converted to 1-based
# row/col, so the river, which sits at Python row (n %/% 2) - 1, is at
# row == n %/% 2 in these tables.
# =============================================================================

suppressPackageStartupMessages({
  library(adegenet); library(dplyr); library(tibble)
})


## ---- Hudson Fst and unbiased heterozygosity for one genlight ---------------
# Returns the full pairwise Fst matrix, per-deme heterozygosity averaged over
# retained SNPs, the number of retained SNPs, and the deme names.
fst_blockwise <- function(gl) {

  pv   <- pop(gl)
  pops <- levels(pv)
  np   <- length(pops)
  nloc <- nLoc(gl)

  # allele frequencies (P) and allele counts (N) per deme per locus
  P <- matrix(NA_real_, np, nloc)
  N <- matrix(NA_real_, np, nloc)
  for (i in seq_len(np)) {
    sub <- as.matrix(gl[which(pv == pops[i]), ])
    nn  <- colSums(!is.na(sub)) * 2       # alleles scored, diploid
    ss  <- colSums(sub, na.rm = TRUE)     # derived allele count
    P[i, ] <- ifelse(nn > 0, ss / nn, NA_real_)
    N[i, ] <- nn
    rm(sub)
  }

  # drop loci that are missing or under-sampled in any deme
  keep <- colSums(is.na(P)) == 0 & colSums(N < 2) == 0
  P <- P[, keep, drop = FALSE]
  N <- N[, keep, drop = FALSE]

  G  <- tcrossprod(P)
  Sq <- rowSums(P^2)
  Tt <- rowSums(P)
  Cc <- rowSums(P * (1 - P) / (N - 1))

  # unbiased expected heterozygosity per deme, averaged over retained SNPs
  hz   <- rowMeans(2 * P * (1 - P) * (N / (N - 1)))
  nloc <- ncol(P)

  num <- outer(Sq, Sq, "+") - 2 * G - outer(Cc, Cc, "+")
  den <- outer(Tt, Tt, "+") - 2 * G
  f   <- num / den
  f[!is.finite(f)] <- NA
  dimnames(f) <- list(pops, pops)

  list(fst = f, hz = setNames(hz, pops), nloc = nloc, pops = pops)
}


## ---- one simulation object, one ecology ------------------------------------
# `x` is a result object from run_both() (it must still carry $meta).
# `sp` is "warm" or "cold".
one_species <- function(x, sp, n = NULL) {

  meta   <- x$meta
  if (is.null(n)) n <- as.integer(meta$params$n)
  seqlen <- as.numeric(meta$params$seq_len)

  r <- fst_blockwise(x[[sp]]$gl_present)

  k   <- as.integer(sub("^p", "", r$pops))
  row <- k %/% n + 1L      # 1-based grid row    (Python i)
  col <- k %%  n + 1L      # 1-based grid column (Python j)

  ij <- which(upper.tri(r$fst), arr.ind = TRUE)

  list(
    pw = tibble(
      run = meta$run, scenario = meta$scenario, ecology = sp,
      deme_i = r$pops[ij[, 1]], deme_j = r$pops[ij[, 2]],
      row_i = row[ij[, 1]], col_i = col[ij[, 1]],
      row_j = row[ij[, 2]], col_j = col[ij[, 2]],
      dist = sqrt((row[ij[, 1]] - row[ij[, 2]])^2 +
                  (col[ij[, 1]] - col[ij[, 2]])^2),
      fst = r$fst[ij],
      grid_n = n),
    dv = tibble(
      run = meta$run, scenario = meta$scenario, ecology = sp,
      deme = r$pops, row = row, col = col,
      # heterozygosity is per retained SNP; spread it over the simulated
      # sequence to get diversity per base pair, comparable with the time series
      pi = as.numeric(r$hz) * r$nloc / seqlen,
      n_snp = r$nloc, seq_len = seqlen,
      grid_n = n))
}


## ---- build the tables over a whole batch -----------------------------------
#
# sims : either a character vector of .rds paths, or a list of already-loaded
#        result objects (e.g. the output of load_batch()).
#
# Genotype matrices are large, so only one ecology is held in memory at a time
# and the garbage collector is called between replicates.
build_tables <- function(sims, verbose = TRUE) {

  PW <- list(); DV <- list()

  for (i in seq_along(sims)) {
    s   <- sims[[i]]
    x   <- if (is.character(s)) readRDS(s) else s
    lbl <- if (is.character(s)) basename(s)
           else if (!is.null(x$meta$file_id)) x$meta$file_id
           else paste0("item", i)

    for (sp in c("cold", "warm")) {
      xs <- list(meta = x$meta); xs[[sp]] <- x[[sp]]
      o  <- one_species(xs, sp)
      PW[[paste(lbl, sp)]] <- o$pw
      DV[[paste(lbl, sp)]] <- o$dv
      rm(xs, o); gc(verbose = FALSE)
    }

    if (is.character(s)) { rm(x); gc(verbose = FALSE) }
    if (verbose) {
      cat(sprintf("  %s done (%d/%d)\n", lbl, i, length(sims)))
      flush.console()
    }
  }

  list(pw = bind_rows(PW), dv = bind_rows(DV))
}


## ---- convenience: run over one or more output directories ------------------
# Writes list(pw = ..., dv = ...) to out_file and returns it.
extract_tables <- function(sim_dirs, out_file = "pw_dv.rds", verbose = TRUE) {

  sims <- unlist(lapply(sim_dirs, list.files,
                        pattern = "^out_.*\\.rds$",
                        full.names = TRUE, recursive = TRUE))
  if (length(sims) == 0)
    stop("No out_*.rds files found in: ", paste(sim_dirs, collapse = ", "))

  cat(sprintf("Extracting Fst and diversity from %d replicate files\n",
              length(sims)))
  res <- build_tables(sims, verbose = verbose)
  saveRDS(res, out_file)
  cat("Written:", out_file, "\n")

  print(res$pw %>% group_by(scenario, ecology) %>%
          summarise(nrun = n_distinct(run), npair = n(),
                    mean_fst = round(mean(fst, na.rm = TRUE), 3),
                    .groups = "drop"))
  print(res$dv %>% group_by(scenario, ecology) %>%
          summarise(nrun = n_distinct(run), ndeme = n(),
                    mean_pi = signif(mean(pi, na.rm = TRUE), 3),
                    .groups = "drop"))

  invisible(res)
}
