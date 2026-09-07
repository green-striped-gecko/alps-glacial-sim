# =============================================================================
# 04_run_parallel.R  --  parallel batch runner
# =============================================================================
#
# Same job as run_batch(), but spread over independent R worker processes.
# Each worker starts its own Python interpreter in the same conda environment
# and runs one design row from start to finish. That is safe because msprime
# carries no state between runs and every row has its own seed and output file.
#
# A PSOCK cluster is used rather than forking: forked processes and an embedded
# Python interpreter do not mix well.
#
# Cost. One replicate (both ecologies, 100 time points) takes roughly one to a
# few hours on one core, so the full 90-row design is a large job. With 30
# workers the published design runs overnight.
#
# This file defines functions only. See scripts/01_run_simulations.R for the
# driver that actually launches a batch.
# =============================================================================

suppressPackageStartupMessages(library(parallel))


## ---- what one worker does ---------------------------------------------------
# Runs in a fresh R process, so it has to set up everything itself: the
# repository path, the Python environment, and the R sources.
.worker_run <- function(row_list, out_dir, home, conda_env) {

  tryCatch({

    Sys.setenv(GLACIALSIM_HOME = home)
    source(file.path(home, "R", "00_setup.R"))
    source(file.path(home, "R", "01_simulate.R"))
    source(file.path(home, "R", "02_design.R"))
    init_python(conda_env, quiet = TRUE)

    rds_path <- file.path(out_dir, paste0(row_list$file_id, ".rds"))
    if (file.exists(rds_path))
      return(list(status = "skipped", file_id = row_list$file_id,
                  file = rds_path))

    t0  <- proc.time()["elapsed"]
    out <- do.call(run_both, design_row_args(row_list))

    out$meta <- list(
      scenario = as.character(row_list$scenario),
      run      = as.integer(row_list$run),
      seed     = as.integer(row_list$seed),
      file_id  = as.character(row_list$file_id),
      params   = row_list[!(names(row_list) %in%
                            c("scenario", "run", "seed", "file_id"))]
    )

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    saveRDS(out, rds_path)

    list(status  = "done",
         file_id = row_list$file_id,
         file    = rds_path,
         elapsed = as.numeric(proc.time()["elapsed"] - t0))

  }, error = function(e)
    list(status = "error", file_id = row_list$file_id,
         message = conditionMessage(e)))
}


## ---- launch a parallel batch ------------------------------------------------
#
# design        design matrix from make_design()
# out_dir       directory for the .rds files
# n_cores       number of worker processes; each one uses a full core and a few
#               hundred MB to a couple of GB of memory, so keep an eye on RAM
# conda_env     conda environment holding msprime; defaults to whatever
#               init_python() would use
# skip_existing rows whose .rds already exists are dropped before launching
run_parallel_batch <- function(design,
                               out_dir       = "batch_output",
                               n_cores       = parallel::detectCores() - 1,
                               conda_env     = NULL,
                               skip_existing = TRUE) {

  home <- glacialsim_home()
  if (is.null(conda_env))
    conda_env <- Sys.getenv("GLACIALSIM_CONDA_ENV", unset = "msprime-env")

  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write.csv(design, file.path(out_dir, "design_matrix.csv"), row.names = FALSE)
  write_session_info(file.path(out_dir, "session_info.txt"))

  if (skip_existing) {
    done <- file.exists(file.path(out_dir, paste0(design$file_id, ".rds")))
    if (any(done)) {
      cat(sprintf("Skipping %d already-completed runs.\n", sum(done)))
      design <- design[!done, , drop = FALSE]
    }
  }

  if (nrow(design) == 0) {
    cat("All runs already complete.\n")
    return(invisible(NULL))
  }

  n_jobs  <- nrow(design)
  n_cores <- min(n_cores, n_jobs)
  cat(sprintf("Launching %d jobs on %d cores...\n", n_jobs, n_cores))

  row_list <- lapply(seq_len(n_jobs), function(i) as.list(design[i, ]))

  cl <- parallel::makeCluster(n_cores, type = "PSOCK",
                              outfile = file.path(out_dir, "worker_log.txt"))
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  # the worker function itself must be shipped to the nodes; everything else it
  # needs is passed as an argument, so no other global state has to be exported
  parallel::clusterExport(cl, ".worker_run", envir = environment(.worker_run))

  t_start <- proc.time()["elapsed"]
  results <- parallel::parLapply(cl, row_list, .worker_run,
                                 out_dir   = out_dir,
                                 home      = home,
                                 conda_env = conda_env)

  statuses <- sapply(results, `[[`, "status")
  cat(sprintf("\n=== Parallel batch finished in %.1f min ===\n",
              (proc.time()["elapsed"] - t_start) / 60))
  cat(sprintf("  done:    %d\n", sum(statuses == "done")))
  cat(sprintf("  skipped: %d\n", sum(statuses == "skipped")))
  cat(sprintf("  errors:  %d\n", sum(statuses == "error")))

  for (e in results[statuses == "error"])
    cat(sprintf("  %s: %s\n", e$file_id, e$message))

  invisible(results)
}
