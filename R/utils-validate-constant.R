validate_project_constant_inputs <- function(cdat, pdat, startp) {
  if (!inherits(cdat, "data.frame") & !inherits(cdat, "matrix")) {
    rlang::abort("\"cdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (!inherits(pdat, "data.frame") & !inherits(pdat, "matrix")) {
    rlang::abort("\"pdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (nrow(cdat) != nrow(pdat)) {
    rlang::abort("\"cdat\" and \"pdat\" must have identical age groups")
  }

  if (ncol(cdat) > ncol(pdat)) {
    rlang::abort(
      "\"pdat\" must include information about all years in \"cdat\""
    )
  } else if (ncol(pdat) == ncol(cdat)) {
    rlang::abort("\"pdat\" must include information in projection years")
  }
  if (ncol(cdat) < 15) {
    rlang::abort("\"cdat\" must have at least 15 years")
  }

  if (ncol(pdat) < 20) {
    rlang::abort("\"pdat\" must have at least 20 years")
  }

  if (!is.null(startp) && !inherits(startp, "numeric")) {
    rlang::abort("\"startp\" must be of type \"numeric\"")
  }

  if (any(pdat[, ncol(cdat)] == 0)) {
    rlang::abort(
      "\"pdat\" has an age group with 0 count in the last known year"
    )
  }
}
