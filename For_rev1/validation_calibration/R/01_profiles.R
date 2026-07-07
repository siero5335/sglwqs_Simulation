load_profiles <- function(path = vc_file("config", "simulation_profiles.yml")) {
  activate_local_lib()
  yaml::read_yaml(path)
}

recursive_modify <- function(x, val) {
  for (nm in names(val)) {
    if (is.list(x[[nm]]) && is.list(val[[nm]])) {
      x[[nm]] <- recursive_modify(x[[nm]], val[[nm]])
    } else {
      x[[nm]] <- val[[nm]]
    }
  }
  x
}

get_profile <- function(profile = Sys.getenv("VALIDATION_CALIBRATION_PROFILE", "smoke")) {
  cfg <- load_profiles()
  if (is.null(cfg$profiles[[profile]])) {
    stop("Unknown validation calibration profile: ", profile, call. = FALSE)
  }
  out <- recursive_modify(cfg$defaults, cfg$profiles[[profile]])
  out$name <- profile
  out$reps <- as.list(out$reps)
  out
}

profile_reps <- function(profile, scenario_family) {
  val <- profile$reps[[scenario_family]]
  if (is.null(val)) {
    stop("Profile does not define reps for scenario family: ", scenario_family, call. = FALSE)
  }
  as.integer(val)
}

profile_scalar <- function(profile, name, type = c("numeric", "integer", "logical", "character")) {
  type <- match.arg(type)
  value <- profile[[name]]
  if (type == "integer") {
    return(as.integer(value))
  }
  if (type == "numeric") {
    return(as.numeric(value))
  }
  if (type == "logical") {
    return(isTRUE(value))
  }
  as.character(value)
}
