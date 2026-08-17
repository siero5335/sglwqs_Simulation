make_base_job <- function(battery,
                          scenario_id,
                          family,
                          n,
                          p,
                          data_seed,
                          method,
                          signal_profile = "baseline",
                          effect_profile = "baseline",
                          group_structure = "legacy13",
                          within_rho = NA_real_,
                          cross_rho = NA_real_,
                          split_replicate = NA_integer_,
                          split_seed = NULL,
                          method_seed = NULL,
                          reuse_from_battery = NA_character_,
                          run_this = TRUE) {
  split_seed <- split_seed %||% (as.integer(data_seed) + 100000L)
  method_seed <- method_seed %||% (as.integer(data_seed) + match(method, r2r2_methods()) * 1000000L)
  data.frame(
    job_id = safe_filename(paste(battery, scenario_id, method, data_seed, split_replicate, sep = "_")),
    battery = battery,
    scenario_id = scenario_id,
    family = family,
    n = as.integer(n),
    p = as.integer(p),
    method = method,
    data_seed = as.integer(data_seed),
    split_seed = as.integer(split_seed),
    method_seed = as.integer(method_seed),
    split_replicate = as.integer(split_replicate),
    signal_profile = signal_profile,
    effect_profile = effect_profile,
    group_structure = group_structure,
    within_rho = as.numeric(within_rho),
    cross_rho = as.numeric(cross_rho),
    reuse_from_battery = reuse_from_battery,
    run_this = isTRUE(run_this),
    stringsAsFactors = FALSE
  )
}

make_sample_size_jobs <- function(seeds = r2r2_seeds(), methods = r2r2_methods()) {
  a1 <- purrr::pmap_dfr(expand.grid(
    n = c(100L, 500L, 1000L, 5000L, 10000L, 50000L),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(n, data_seed, method) {
    make_base_job(
      battery = "A1_gaussian_sample_size",
      scenario_id = sprintf("A1_gaussian_n%s_baseline", n),
      family = "gaussian",
      n = n,
      p = 13L,
      data_seed = data_seed,
      method = method,
      signal_profile = "baseline",
      effect_profile = "baseline",
      group_structure = "legacy13"
    )
  })
  a2 <- purrr::pmap_dfr(expand.grid(
    n = c(500L, 1000L, 5000L, 10000L),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(n, data_seed, method) {
    make_base_job(
      battery = "A2_weak_gaussian_sample_size",
      scenario_id = sprintf("A2_gaussian_n%s_weak", n),
      family = "gaussian",
      n = n,
      p = 13L,
      data_seed = data_seed,
      method = method,
      signal_profile = "weak",
      effect_profile = "weak",
      group_structure = "legacy13"
    )
  })
  dplyr::bind_rows(a1, a2)
}

make_highdim_jobs <- function(seeds = r2r2_seeds(), methods = r2r2_methods()) {
  b1 <- purrr::pmap_dfr(expand.grid(
    p = c(50L, 100L, 200L),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(p, data_seed, method) {
    make_base_job(
      battery = "B1_gaussian_highdim",
      scenario_id = sprintf("B1_gaussian_n1000_p%03d_baseline", p),
      family = "gaussian",
      n = 1000L,
      p = p,
      data_seed = data_seed,
      method = method,
      signal_profile = "baseline",
      effect_profile = "baseline",
      group_structure = "equal10",
      within_rho = 0.45,
      cross_rho = 0.05
    )
  })
  b2 <- purrr::pmap_dfr(expand.grid(
    p = c(50L, 100L, 200L),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), function(p, data_seed, method) {
    make_base_job(
      battery = "B2_weak_gaussian_highdim",
      scenario_id = sprintf("B2_gaussian_n1000_p%03d_weak", p),
      family = "gaussian",
      n = 1000L,
      p = p,
      data_seed = data_seed,
      method = method,
      signal_profile = "weak",
      effect_profile = "weak",
      group_structure = "equal10",
      within_rho = 0.45,
      cross_rho = 0.05
    )
  })
  dplyr::bind_rows(b1, b2)
}

make_split_stability_jobs <- function(split_reps = 100L) {
  scenarios <- data.frame(
    scenario_label = c("global_null_n500", "global_null_n5000", "partial_null_n500", "partial_null_n5000"),
    effect_profile = c("global_null", "global_null", "partial_null", "partial_null"),
    n = c(500L, 5000L, 500L, 5000L),
    within_rho = c(0.20, 0.20, 0.60, 0.60),
    cross_rho = c(0, 0, 0.30, 0.30),
    data_seed = c(991001L, 991002L, 991101L, 991102L),
    stringsAsFactors = FALSE
  )
  purrr::pmap_dfr(scenarios, function(scenario_label, effect_profile, n, within_rho, cross_rho, data_seed) {
    purrr::map_dfr(seq_len(split_reps), function(split_id) {
      make_base_job(
        battery = "C_gaussian_split_stability",
        scenario_id = sprintf("C_gaussian_split_%s", scenario_label),
        family = "gaussian",
        n = n,
        p = 13L,
        data_seed = data_seed,
        method = "SGL-WQS",
        signal_profile = if (identical(effect_profile, "global_null")) "global_null" else "baseline",
        effect_profile = effect_profile,
        group_structure = "legacy13",
        within_rho = within_rho,
        cross_rho = cross_rho,
        split_replicate = split_id,
        split_seed = as.integer(data_seed) + 70000L + split_id,
        method_seed = as.integer(data_seed) + 90000L + split_id
      )
    })
  })
}

make_factorial_jobs <- function(seeds = r2r2_seeds(), methods = r2r2_methods()) {
  grid <- expand.grid(
    family = c("binomial", "gaussian"),
    group_structure = c("balanced", "unbalanced"),
    effect_profile = c("uniform_baseline", "heterogeneous", "weak_heterogeneous"),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  purrr::pmap_dfr(grid, function(family, group_structure, effect_profile, data_seed, method) {
    make_base_job(
      battery = "D_unbalanced_effect_factorial",
      scenario_id = sprintf("D_%s_%s_%s_n1000_p100", family, group_structure, effect_profile),
      family = family,
      n = 1000L,
      p = 100L,
      data_seed = data_seed,
      method = method,
      signal_profile = effect_profile,
      effect_profile = effect_profile,
      group_structure = group_structure,
      within_rho = 0.45,
      cross_rho = 0.05
    )
  })
}

make_hard_setting_jobs <- function(seeds = r2r2_seeds(), methods = r2r2_methods()) {
  run_grid <- expand.grid(
    family = c("binomial", "gaussian"),
    n = c(500L, 5000L),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  run_jobs <- purrr::pmap_dfr(run_grid, function(family, n, data_seed, method) {
    make_base_job(
      battery = "E_hard_setting_sample_size",
      scenario_id = sprintf("E_%s_unbalanced_weak_heterogeneous_n%s_p100", family, n),
      family = family,
      n = n,
      p = 100L,
      data_seed = data_seed,
      method = method,
      signal_profile = "weak_heterogeneous",
      effect_profile = "weak_heterogeneous",
      group_structure = "unbalanced",
      within_rho = 0.45,
      cross_rho = 0.05
    )
  })

  reuse_grid <- expand.grid(
    family = c("binomial", "gaussian"),
    data_seed = seeds,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  reuse_jobs <- purrr::pmap_dfr(reuse_grid, function(family, data_seed, method) {
    make_base_job(
      battery = "E_hard_setting_sample_size",
      scenario_id = sprintf("E_%s_unbalanced_weak_heterogeneous_n1000_p100", family),
      family = family,
      n = 1000L,
      p = 100L,
      data_seed = data_seed,
      method = method,
      signal_profile = "weak_heterogeneous",
      effect_profile = "weak_heterogeneous",
      group_structure = "unbalanced",
      within_rho = 0.45,
      cross_rho = 0.05,
      reuse_from_battery = "D_unbalanced_effect_factorial",
      run_this = FALSE
    )
  })
  dplyr::bind_rows(run_jobs, reuse_jobs)
}

make_production_manifest <- function() {
  dplyr::bind_rows(
    make_sample_size_jobs(),
    make_highdim_jobs(),
    make_split_stability_jobs(100L),
    make_factorial_jobs(),
    make_hard_setting_jobs()
  ) |>
    dplyr::arrange(.data$battery, .data$scenario_id, .data$data_seed, .data$method, .data$split_replicate)
}

make_smoke_manifest <- function() {
  smoke_seeds <- r2r2_seeds(2L)
  methods <- r2r2_methods()
  dplyr::bind_rows(
    make_base_job("A1_gaussian_sample_size", "smoke_A_gaussian_n100", "gaussian", 100L, 13L, smoke_seeds[1], methods),
    make_base_job("A1_gaussian_sample_size", "smoke_A_binomial_n500", "binomial", 500L, 13L, smoke_seeds[2], methods),
    make_base_job("B1_gaussian_highdim", "smoke_B_gaussian_p200", "gaussian", 1000L, 200L, smoke_seeds[1], methods,
      signal_profile = "baseline", effect_profile = "baseline", group_structure = "equal10", within_rho = 0.45, cross_rho = 0.05
    ),
    make_base_job("C_gaussian_split_stability", "smoke_C_gaussian_split_global_n500", "gaussian", 500L, 13L, 991001L, "SGL-WQS",
      signal_profile = "global_null", effect_profile = "global_null", group_structure = "legacy13",
      within_rho = 0.20, cross_rho = 0, split_replicate = 1L, split_seed = 991001L + 70001L, method_seed = 991001L + 90001L
    ),
    purrr::pmap_dfr(expand.grid(
      family = c("binomial", "gaussian"),
      group_structure = c("balanced", "unbalanced"),
      method = c("SGL-WQS", "qgcomp"),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ), function(family, group_structure, method) {
      make_base_job(
        "D_unbalanced_effect_factorial",
        sprintf("smoke_D_%s_%s_heterogeneous", family, group_structure),
        family, 500L, 100L, smoke_seeds[2], method,
        signal_profile = "heterogeneous",
        effect_profile = "heterogeneous",
        group_structure = group_structure,
        within_rho = 0.45,
        cross_rho = 0.05
      )
    }),
    purrr::pmap_dfr(data.frame(
      family = c("binomial", "gaussian", "gaussian"),
      method = c("qgcomp", "SGL-WQS", "qgcomp"),
      stringsAsFactors = FALSE
    ), function(family, method) {
      make_base_job(
        "E_hard_setting_sample_size",
        sprintf("smoke_E_%s_unbalanced_weak_n500", family),
        family, 500L, 100L, smoke_seeds[1], method,
        signal_profile = "weak_heterogeneous",
        effect_profile = "weak_heterogeneous",
        group_structure = "unbalanced",
        within_rho = 0.45,
        cross_rho = 0.05
      )
    })
  ) |>
    dplyr::mutate(job_id = paste0("smoke_", .data$job_id))
}

write_manifest <- function(manifest, name = "scenario_manifest.csv") {
  r2r2_make_dirs()
  atomic_write_csv(manifest, r2r2_file("config", name))
}

manifest_job_paths <- function(manifest) {
  vapply(split(manifest, seq_len(nrow(manifest))), job_output_path, character(1))
}

run_job_table <- function(jobs,
                          settings = r2r2_settings(),
                          workers = 1L,
                          force = FALSE) {
  jobs <- jobs[jobs$run_this, , drop = FALSE]
  if (!nrow(jobs)) {
    return(list())
  }
  jobs$path <- manifest_job_paths(jobs)
  todo <- jobs[force | !vapply(jobs$path, is_complete_job_file, logical(1)), , drop = FALSE]
  message("Selected jobs: ", nrow(jobs), "; already complete: ", nrow(jobs) - nrow(todo), "; to run: ", nrow(todo))
  if (!nrow(todo)) {
    return(atomic_read_results(jobs$path))
  }
  job_list <- split(todo, seq_len(nrow(todo)))
  worker_fun <- function(j) {
    source(file.path(Sys.getenv("R2R2_ROOT", unset = file.path(getwd(), "reviewer2_round2")), "R", "design_helpers.R"), chdir = TRUE)
    source_r2r2("gaussian_matched_generators.R")
    source_r2r2("unbalanced_factorial_generator.R")
    source_r2r2("method_wrappers.R")
    r2r2_set_thread_env()
    run_method_job_resumable(j, settings = settings, force = force)
  }

  if (as.integer(workers) > 1L && .Platform$OS.type == "unix") {
    par_res <- tryCatch(
      parallel::mclapply(job_list, worker_fun, mc.cores = as.integer(workers), mc.preschedule = FALSE),
      error = function(e) e
    )
    if (inherits(par_res, "error")) {
      message("Parallel backend failed; retrying the same jobs sequentially: ", conditionMessage(par_res))
      par_res <- lapply(job_list, worker_fun)
    } else {
      failed_parallel <- vapply(par_res, inherits, logical(1), what = "try-error")
      if (any(failed_parallel)) {
        message("Some forked jobs returned try-error; retrying those jobs sequentially.")
        retry <- lapply(job_list[failed_parallel], worker_fun)
        par_res[failed_parallel] <- retry
      }
    }
  } else {
    par_res <- lapply(job_list, worker_fun)
  }
  invisible(par_res)
  atomic_read_results(jobs$path)
}

summarize_manifest_status <- function(manifest) {
  run_jobs <- manifest[manifest$run_this, , drop = FALSE]
  paths <- manifest_job_paths(run_jobs)
  status <- run_jobs |>
    dplyr::mutate(
      output_path = paths,
      completed = file.exists(paths),
      valid_complete = vapply(paths, is_complete_job_file, logical(1)),
      size_bytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_)
    )
  atomic_write_csv(status, r2r2_result_file("summaries", "job_completion_status.csv"))
  status
}
