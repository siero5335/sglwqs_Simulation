dense_group_signal_batteries <- function(include_smoke = FALSE) {
  production <- "H_dense_group_signal"
  if (include_smoke) c(production, "SMOKE_H_dense_group_signal") else production
}

is_dense_group_signal_battery <- function(battery) {
  battery %in% dense_group_signal_batteries(include_smoke = TRUE)
}

dense_group_signal_selected_families <- function() {
  raw <- Sys.getenv(
    "R2R2_DENSE_FAMILIES",
    unset = Sys.getenv("R2R2_DENSE_FAMILY", unset = "binomial,gaussian")
  )
  families <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  families <- unique(families[nzchar(families)])
  if (identical(families, "all")) families <- c("binomial", "gaussian")
  unknown <- setdiff(families, c("binomial", "gaussian"))
  if (!length(families) || length(unknown)) {
    stop(
      "R2R2_DENSE_FAMILIES must contain binomial, gaussian, or all. Unknown: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  families
}

make_dense_group_signal_manifest <- function(seeds = r2r2_seeds(),
                                             methods = r2r2_methods(),
                                             sample_sizes = c(500L, 5000L)) {
  grid <- expand.grid(
    family = c("binomial", "gaussian"),
    n = sample_sizes,
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  purrr::pmap_dfr(grid, function(family, n, data_seed, method) {
    make_base_job(
      battery = "H_dense_group_signal",
      scenario_id = sprintf("H_%s_unbalanced_dense_group_n%s_p100", family, n),
      family = family,
      n = n,
      p = 100L,
      data_seed = data_seed,
      method = method,
      signal_profile = "dense_group_matched",
      effect_profile = "dense_group_matched",
      group_structure = "unbalanced",
      within_rho = 0.45,
      cross_rho = 0.05
    )
  }) |>
    dplyr::arrange(.data$family, .data$n, .data$data_seed, .data$method)
}

make_dense_group_signal_smoke_manifest <- function() {
  seed <- r2r2_seeds(1L)[[1L]]
  purrr::pmap_dfr(expand.grid(
    family = c("binomial", "gaussian"),
    method = r2r2_methods(),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(family, method) {
    make_base_job(
      battery = "SMOKE_H_dense_group_signal",
      scenario_id = sprintf("smoke_H_%s_unbalanced_dense_group_n500_p100", family),
      family = family,
      n = 500L,
      p = 100L,
      data_seed = seed,
      method = method,
      signal_profile = "dense_group_matched",
      effect_profile = "dense_group_matched",
      group_structure = "unbalanced",
      within_rho = 0.45,
      cross_rho = 0.05
    )
  })
}

generate_dense_group_signal_job_data <- function(job) {
  generate_dense_group_signal_data(
    n = as.integer(job$n),
    data_seed = as.integer(job$data_seed),
    family = as.character(job$family),
    within_rho = as.numeric(job$within_rho),
    cross_rho = as.numeric(job$cross_rho),
    gaussian_sd = 1
  )
}

install_dense_group_signal_dispatch <- function(env = .GlobalEnv) {
  current <- get("generate_job_data", envir = env, inherits = FALSE)
  if (isTRUE(attr(current, "dense_group_signal_dispatch"))) {
    return(invisible(current))
  }
  base_generator <- current
  dispatch <- function(job) {
    if (is_dense_group_signal_battery(as.character(job$battery))) {
      return(generate_dense_group_signal_job_data(job))
    }
    base_generator(job)
  }
  attr(dispatch, "dense_group_signal_dispatch") <- TRUE
  assign("generate_job_data", dispatch, envir = env)
  invisible(dispatch)
}

dense_group_signal_manifest_paths <- function(manifest) {
  vapply(split(manifest, seq_len(nrow(manifest))), job_output_path, character(1))
}

run_dense_group_signal_jobs <- function(jobs,
                                        settings = r2r2_settings(),
                                        workers = 1L,
                                        force = FALSE) {
  jobs <- jobs[jobs$run_this, , drop = FALSE]
  jobs$path <- dense_group_signal_manifest_paths(jobs)
  todo <- jobs[force | !vapply(jobs$path, is_complete_job_file, logical(1)), , drop = FALSE]
  message(
    "Selected dense group-signal jobs: ", nrow(jobs),
    "; already complete: ", nrow(jobs) - nrow(todo),
    "; to run: ", nrow(todo)
  )
  if (!nrow(todo)) return(atomic_read_results(jobs$path))

  job_list <- split(todo, seq_len(nrow(todo)))
  worker_fun <- function(j) {
    source(file.path(
      Sys.getenv("R2R2_ROOT", unset = file.path(getwd(), "reviewer2_round2")),
      "R", "design_helpers.R"
    ), chdir = TRUE)
    source_r2r2("gaussian_matched_generators.R")
    source_r2r2("unbalanced_factorial_generator.R")
    source_r2r2("method_wrappers.R")
    source_r2r2("dense_group_signal_generator.R")
    source_r2r2("dense_group_signal_helpers.R")
    install_dense_group_signal_dispatch(.GlobalEnv)
    r2r2_set_thread_env()
    run_method_job_resumable(j, settings = settings, force = force)
  }

  if (as.integer(workers) > 1L && .Platform$OS.type == "unix") {
    results <- tryCatch(
      parallel::mclapply(job_list, worker_fun,
        mc.cores = as.integer(workers), mc.preschedule = FALSE
      ),
      error = function(e) e
    )
    if (inherits(results, "error")) {
      message("Parallel backend failed; retrying all selected jobs sequentially: ", conditionMessage(results))
      results <- lapply(job_list, worker_fun)
    } else {
      failed_parallel <- vapply(results, inherits, logical(1), what = "try-error")
      if (any(failed_parallel)) {
        message("Forked jobs returned try-error; retrying those jobs sequentially.")
        results[failed_parallel] <- lapply(job_list[failed_parallel], worker_fun)
      }
    }
  } else {
    results <- lapply(job_list, worker_fun)
  }
  invisible(results)
  atomic_read_results(jobs$path)
}

summarize_dense_group_signal_status <- function(manifest,
                                                filename = "dense_group_signal_completion_status.csv") {
  paths <- dense_group_signal_manifest_paths(manifest)
  status <- manifest |>
    dplyr::mutate(
      output_path = paths,
      completed = file.exists(paths),
      valid_complete = vapply(paths, is_complete_job_file, logical(1)),
      size_bytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_)
    )
  atomic_write_csv(status, r2r2_result_file("summaries", filename))
  status
}
