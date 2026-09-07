# =============================================================================
# 03_run_batch.R  --  serial batch runner and result loaders
# =============================================================================
#
# run_batch() works through a design matrix one row at a time and writes one
# .rds file per row. It is the simple, dependency-free way to run the study;
# use 04_run_parallel.R on a multi-core machine.
#
# Output layout
#   <out_dir>/design_matrix.csv     the design that was run
#   <out_dir>/session_info.txt      R and Python versions used
#   <out_dir>/out_<scenario>_run<RR>.rds
#
# Each .rds holds the list returned by run_both(), i.e. $warm and $cold, plus
# a $meta entry recording the scenario, replicate number, seed and parameters.
# =============================================================================


## ---- run every row of a design ---------------------------------------------
#
# design        design matrix from make_design()
# out_dir       directory for the .rds files, created if needed
# skip_existing TRUE lets an interrupted batch be resumed safely
run_batch <- function(design,
                      out_dir       = "batch_output",
                      skip_existing = TRUE) {

  if (!exists("gs"))
    stop("Python not initialised - call init_python() first (see R/00_setup.R)")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat("Created output directory:", out_dir, "\n")
  }

  write.csv(design, file.path(out_dir, "design_matrix.csv"), row.names = FALSE)
  write_session_info(file.path(out_dir, "session_info.txt"))

  n_total <- nrow(design)
  t_start <- proc.time()["elapsed"]

  for (i in seq_len(n_total)) {
    row      <- design[i, ]
    rds_path <- file.path(out_dir, paste0(row$file_id, ".rds"))

    cat(sprintf("\n[%d/%d] scenario=%s  run=%02d  seed=%d\n",
                i, n_total, row$scenario, row$run, row$seed))

    if (skip_existing && file.exists(rds_path)) {
      cat("  -> file exists, skipping.\n")
      next
    }

    t0  <- proc.time()["elapsed"]
    out <- tryCatch(do.call(run_both, design_row_args(row)),
                    error = function(e) {
                      message("  ERROR in run_both(): ", conditionMessage(e))
                      NULL
                    })

    if (is.null(out)) {
      cat("  -> run failed, nothing saved.\n")
      next
    }

    out$meta <- list(
      scenario = as.character(row$scenario),
      run      = as.integer(row$run),
      seed     = as.integer(row$seed),
      file_id  = as.character(row$file_id),
      params   = as.list(row[, !(names(row) %in%
                                 c("scenario", "run", "seed", "file_id"))])
    )

    saveRDS(out, rds_path)

    elapsed <- proc.time()["elapsed"] - t0
    eta     <- (proc.time()["elapsed"] - t_start) / i * (n_total - i)
    cat(sprintf("  -> saved %s  (%.1f min; about %.0f min left)\n",
                rds_path, elapsed / 60, eta / 60))
  }

  cat("\n=== Batch complete ===\n")
  invisible(design)
}


## ---- load a finished batch --------------------------------------------------
#
# Returns a named list of result objects, named by file_id, e.g.
#   results[["out_1peak_run03"]]$cold$pi
#
# out_dir   : one directory, or several - they are read in the order given
# scenarios : NULL for everything, or a subset such as c("1peak", "3peak")
load_batch <- function(out_dir = "batch_output", scenarios = NULL) {

  rds_files <- unlist(lapply(out_dir, list.files,
                             pattern = "^out_.*\\.rds$",
                             full.names = TRUE, recursive = TRUE))

  if (!is.null(scenarios)) {
    pat <- sprintf("out_(%s)_run", paste(scenarios, collapse = "|"))
    rds_files <- rds_files[grepl(pat, basename(rds_files))]
  }

  if (length(rds_files) == 0) {
    warning("No matching .rds files found in: ", paste(out_dir, collapse = ", "))
    return(list())
  }

  cat(sprintf("Loading %d replicate files ...\n", length(rds_files)))
  results <- lapply(rds_files, readRDS)
  names(results) <- tools::file_path_sans_ext(basename(rds_files))
  cat("Done.\n")
  results
}


## ---- pull one statistic out of a loaded batch ------------------------------
#
# Returns a matrix of time points x replicates.
#
# scenario : "1peak" | "3peak" | "3peak_river"
# species  : "warm"  | "cold"
# stat     : "pi" | "he_mean" | "fst_mean" | "tajima_d" | "theta_w"
extract_stat_matrix <- function(results, scenario, species, stat = "pi") {
  keep <- grepl(paste0("^out_", scenario, "_run"), names(results))
  reps <- results[keep]
  if (length(reps) == 0) stop("No results found for scenario: ", scenario)
  mat <- do.call(cbind, lapply(reps, function(r) r[[species]]$pi[[stat]]))
  colnames(mat) <- names(reps)
  mat
}
