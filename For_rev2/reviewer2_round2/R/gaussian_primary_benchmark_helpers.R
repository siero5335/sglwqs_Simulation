gaussian_primary_benchmark_batteries <- function(include_smoke = FALSE) {
  production <- c(
    "F_gaussian_correlation_robustness",
    "G_gaussian_active_component"
  )
  if (!include_smoke) {
    return(production)
  }
  c(
    production,
    "SMOKE_F_gaussian_correlation_robustness",
    "SMOKE_G_gaussian_active_component"
  )
}

is_gaussian_primary_benchmark_battery <- function(battery) {
  battery %in% gaussian_primary_benchmark_batteries(include_smoke = TRUE)
}

make_gaussian_primary_benchmark_manifest <- function(seeds = r2r2_seeds(),
                                                     methods = r2r2_methods(),
                                                     n = 10000L,
                                                     correlations = c(0.2, 0.5, 0.8, 0.95)) {
  correlation_jobs <- purrr::pmap_dfr(expand.grid(
    correlation = correlations,
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(correlation, data_seed, method) {
    make_base_job(
      battery = "F_gaussian_correlation_robustness",
      scenario_id = sprintf("F_gaussian_correlation_r%s_n%s", format(correlation, trim = TRUE), n),
      family = "gaussian",
      n = n,
      p = 13L,
      data_seed = data_seed,
      method = method,
      signal_profile = "sparse_correlation",
      effect_profile = "sparse_correlation",
      group_structure = "legacy13_correlation",
      within_rho = correlation,
      cross_rho = 0
    )
  })

  active_jobs <- purrr::pmap_dfr(expand.grid(
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(data_seed, method) {
    make_base_job(
      battery = "G_gaussian_active_component",
      scenario_id = sprintf("G_gaussian_single_active_component_n%s", n),
      family = "gaussian",
      n = n,
      p = 13L,
      data_seed = data_seed,
      method = method,
      signal_profile = "single_active_component",
      effect_profile = "single_active_component",
      group_structure = "legacy13_active_component",
      within_rho = 0.70,
      cross_rho = 0
    )
  })

  dplyr::bind_rows(correlation_jobs, active_jobs) |>
    dplyr::arrange(.data$battery, .data$within_rho, .data$data_seed, .data$method)
}

make_gaussian_primary_benchmark_smoke_manifest <- function() {
  methods <- r2r2_methods()
  seed <- r2r2_seeds(1L)[[1L]]
  correlation_jobs <- purrr::pmap_dfr(expand.grid(
    correlation = c(0.2, 0.95),
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(correlation, method) {
    make_base_job(
      battery = "SMOKE_F_gaussian_correlation_robustness",
      scenario_id = sprintf("smoke_F_gaussian_correlation_r%s_n500", format(correlation, trim = TRUE)),
      family = "gaussian",
      n = 500L,
      p = 13L,
      data_seed = seed,
      method = method,
      signal_profile = "sparse_correlation",
      effect_profile = "sparse_correlation",
      group_structure = "legacy13_correlation",
      within_rho = correlation,
      cross_rho = 0
    )
  })
  active_jobs <- purrr::map_dfr(methods, function(method) {
    make_base_job(
      battery = "SMOKE_G_gaussian_active_component",
      scenario_id = "smoke_G_gaussian_single_active_component_n500",
      family = "gaussian",
      n = 500L,
      p = 13L,
      data_seed = seed,
      method = method,
      signal_profile = "single_active_component",
      effect_profile = "single_active_component",
      group_structure = "legacy13_active_component",
      within_rho = 0.70,
      cross_rho = 0
    )
  })
  dplyr::bind_rows(correlation_jobs, active_jobs)
}

generate_gaussian_primary_benchmark_job_data <- function(job) {
  if (!identical(as.character(job$family), "gaussian")) {
    stop("Gaussian primary benchmark jobs must use family='gaussian'.", call. = FALSE)
  }
  if (grepl("correlation_robustness", as.character(job$battery), fixed = TRUE)) {
    return(generate_gaussian_correlation_benchmark_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      correlation = as.numeric(job$within_rho),
      gaussian_sd = 1
    ))
  }
  if (grepl("active_component", as.character(job$battery), fixed = TRUE)) {
    return(generate_gaussian_active_component_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      pcb_correlation = as.numeric(job$within_rho),
      gaussian_sd = 1
    ))
  }
  stop("Unknown Gaussian primary benchmark battery: ", job$battery, call. = FALSE)
}

install_gaussian_primary_benchmark_dispatch <- function(env = .GlobalEnv) {
  current <- get("generate_job_data", envir = env, inherits = FALSE)
  if (isTRUE(attr(current, "gaussian_primary_benchmark_dispatch"))) {
    return(invisible(current))
  }
  base_generator <- current
  dispatch <- function(job) {
    if (is_gaussian_primary_benchmark_battery(as.character(job$battery))) {
      return(generate_gaussian_primary_benchmark_job_data(job))
    }
    base_generator(job)
  }
  attr(dispatch, "gaussian_primary_benchmark_dispatch") <- TRUE
  assign("generate_job_data", dispatch, envir = env)
  invisible(dispatch)
}

gaussian_primary_manifest_paths <- function(manifest) {
  vapply(split(manifest, seq_len(nrow(manifest))), job_output_path, character(1))
}

run_gaussian_primary_benchmark_jobs <- function(jobs,
                                                settings = r2r2_settings(),
                                                workers = 1L,
                                                force = FALSE) {
  jobs <- jobs[jobs$run_this, , drop = FALSE]
  if (!nrow(jobs)) {
    return(list())
  }
  jobs$path <- gaussian_primary_manifest_paths(jobs)
  todo <- jobs[force | !vapply(jobs$path, is_complete_job_file, logical(1)), , drop = FALSE]
  message(
    "Selected Gaussian primary benchmark jobs: ", nrow(jobs),
    "; already complete: ", nrow(jobs) - nrow(todo),
    "; to run: ", nrow(todo)
  )
  if (!nrow(todo)) {
    return(atomic_read_results(jobs$path))
  }

  job_list <- split(todo, seq_len(nrow(todo)))
  worker_fun <- function(j) {
    source(file.path(
      Sys.getenv("R2R2_ROOT", unset = file.path(getwd(), "reviewer2_round2")),
      "R", "design_helpers.R"
    ), chdir = TRUE)
    source_r2r2("gaussian_matched_generators.R")
    source_r2r2("unbalanced_factorial_generator.R")
    source_r2r2("method_wrappers.R")
    source_r2r2("gaussian_primary_benchmark_generators.R")
    source_r2r2("gaussian_primary_benchmark_helpers.R")
    install_gaussian_primary_benchmark_dispatch(.GlobalEnv)
    r2r2_set_thread_env()
    run_method_job_resumable(j, settings = settings, force = force)
  }

  if (as.integer(workers) > 1L && .Platform$OS.type == "unix") {
    results <- tryCatch(
      parallel::mclapply(
        job_list,
        worker_fun,
        mc.cores = as.integer(workers),
        mc.preschedule = FALSE
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

summarize_gaussian_primary_benchmark_status <- function(manifest,
                                                        filename = "gaussian_primary_benchmark_completion_status.csv") {
  paths <- gaussian_primary_manifest_paths(manifest)
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
