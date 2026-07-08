#' Validate inputs for summary functions
#'
#' @keywords internal
validate_summary_inputs <- function(
  printpred = NULL,
  printcall = NULL,
  digits = NULL
) {
  if (!is.null(printpred) && !inherits(printpred, "logical")) {
    stop("Variable \"printpred\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!is.null(printcall) && !inherits(printcall, "logical")) {
    stop("Variable \"printcall\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!is.null(digits)) {
    if (!inherits(digits, "numeric")) {
      stop("Variable \"digits\" must be of type \"numeric\"")
    }

    if (digits < 0) {
      stop("Variable \"digits\" must be >= 0")
    }
  }
}

#' Validate inputs for getproj functions
#'
#' @keywords internal
validate_getproj_inputs <- function(
  object,
  standpop = NULL,
  pdat = NULL,
  cdat = NULL,
  startp = NULL
) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)

    if (nrow(pdat) != length(standpop@weights)) {
      stop(
        "Variable \"standpop\" must have same number of age groups as \"pdat\""
      )
    }
  }

  if (!inherits(pdat, "data.frame") & !inherits(pdat, "matrix")) {
    stop("\"pdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (!is.null(object$predictions)) {
    if (nrow(pdat) != nrow(object$predictions)) {
      stop(
        "\"pdat\" must have the same number of age groups as predictions"
      )
    }
  } else {
    if (nrow(pdat) != nrow(object$agsproj)) {
      stop(
        "\"pdat\" must have same number of age groups as predictions"
      )
    }
  }

  if (!is.null(cdat)) {
    if (!inherits(cdat, "data.frame") & !inherits(cdat, "matrix")) {
      stop("\"cdat\" must be of type \"data.frame\" or \"matrix\"")
    }

    if (nrow(cdat) != nrow(object$predictions)) {
      stop(
        "\"cdat\" must have the same number of age groups as predictions"
      )
    }
  }

  if (!is.null(startp)) {
    if (!inherits(startp, "numeric")) {
      stop("\"startp\" must be of type \"numeric\"")
    }
  }
}

#' Validate inputs for getpred functions
#'
#' @keywords internal
validate_getpred_inputs <- function(
  byage,
  standpop,
  incidence,
  agegroups,
  excludeobs
) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  if (!inherits(incidence, "logical")) {
    stop("Variable \"incidence\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!inherits(excludeobs, "logical")) {
    stop("Variable \"excludeobs\" must be \"TRUE\" or \"FALSE\"")
  }

  if ((!is.null(standpop)) && (!incidence)) {
    stop(
      "\"standpop\" should only be used with incidence predictions (incidence=T)"
    )
  }

  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)

    if ((length(standpop) != length(agegroups)) && (agegroups[1] != "all")) {
      stop("\"standpop\" must be the same length as \"agegroups\"")
    }

    if (byage) {
      stop("\"standpop\" is only valid for \"byage=F\"")
    }
  }
}

#' Validate inputs for prediction functions
#'
#' @keywords internal
validate_prediction_inputs <- function(cuttrd, shortp) {
  if (!inherits(cuttrd, "numeric")) {
    stop("Variable \"cuttrd\" must be of type \"numeric\"")
  }

  if (cuttrd > 1 | cuttrd < 0) {
    stop("Variable \"cuttrd\" must be >= 0 and <= 1")
  }

  if (!inherits(shortp, "numeric")) {
    stop("Variable \"shortp\" must be of type \"numeric\"")
  }

  if (shortp > 1 | shortp < 0) {
    stop("Variable \"shortp\" must be >= 0 and <= 1")
  }
}

#' Validate inputs for estimate functions
#'
#' @keywords internal
validate_estimate_inputs <- function(pyr, cases, pGOF, linkfunc) {
  if (ncol(cases) > ncol(pyr)) {
    stop("\"pyr\" must include information about all periods in \"cases\"")
  }

  if (ncol(pyr) == ncol(cases)) {
    stop("\"pyr\" must include information on future rates")
  }

  if (nrow(cases) != nrow(pyr)) {
    stop("\"cases\" and \"pyr\" must have the same number of age groups")
  }

  if (ncol(cases) < 2) {
    stop("\"cases\" must have at least 2 age groups")
  }

  if (ncol(pyr) < 2) {
    stop("\"pyr\" must have at least 2 age groups")
  }

  if (!inherits(pGOF, "numeric")) {
    stop("Variable \"pGOF\" must be of type \"numeric\"")
  }

  if (pGOF > 1 | pGOF < 0) {
    stop("Variable \"pGOF\" must be >= 0 and <= 1")
  }

  if (!(linkfunc) %in% list("power5", "sqrt", "log", "identity")) {
    stop(
      "Variable \"linkfunc\" must be one of \"power5\", \"log\", \"sqrt\", or \"identity\""
    )
  }
}

#' Validate inputs for projection functions
#'
#' @keywords internal
validate_projection_inputs <- function(
  cdat,
  pdat,
  projfor = NULL,
  n5case = NULL,
  linkfunc = NULL,
  startp = NULL,
  newcohort = NULL,
  pGOF = NULL,
  pD = NULL,
  shortp = NULL,
  cuttrd = NULL,
  nagg = NULL,
  startage = NULL,
  ncase = NULL
) {
  if (!inherits(cdat, "data.frame") & !inherits(cdat, "matrix")) {
    stop("\"cdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (!inherits(pdat, "data.frame") & !inherits(pdat, "matrix")) {
    stop("\"pdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (nrow(cdat) != nrow(pdat)) {
    stop("\"cdat\" and \"pdat\" must have identical age groups")
  }

  if (ncol(cdat) > ncol(pdat)) {
    stop("\"pdat\" must include information about all years in \"cdat\"")
  } else if (ncol(pdat) == ncol(cdat)) {
    stop("\"pdat\" must include information in projection years")
  }

  if (!is.null(n5case) && !inherits(n5case, "numeric")) {
    stop("\"n5case\" must be of type \"numeric\" or NULL")
  }

  if (!is.null(projfor) && !(projfor %in% list("incidence", "mortality"))) {
    stop("\"projfor\" must be one of \"incidence\" or \"mortality\"")
  }

  if (
    !is.null(linkfunc) &&
      !(linkfunc) %in% list("power5", "sqrt", "log", "identity")
  ) {
    stop(
      "Variable \"linkfunc\" must be one of \"power5\", \"log\", \"sqrt\", or \"identity\""
    )
  }

  if (!is.null(startp) && !inherits(startp, "numeric")) {
    stop("\"startp\" must be of type \"numeric\"")
  }

  if (!is.null(newcohort) && !inherits(newcohort, "logical")) {
    stop("\"newcohort\" must be of type \"logical\"")
  }

  if (!is.null(pGOF)) {
    if (!inherits(pGOF, "numeric")) {
      stop("Variable \"pGOF\" must be of type \"numeric\"")
    }

    if (pGOF > 1 | pGOF < 0) {
      stop("Variable \"pGOF\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(pD)) {
    if (!inherits(pD, "numeric")) {
      stop("Variable \"pD\" must be of type \"numeric\"")
    }

    if (pD > 1 | pD < 0) {
      stop("Variable \"pD\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(shortp)) {
    if (!inherits(shortp, "numeric")) {
      stop("Variable \"shortp\" must be of type \"numeric\"")
    }

    if (shortp > 1 | shortp < 0) {
      stop("Variable \"shortp\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(cuttrd)) {
    if (!inherits(cuttrd, "numeric")) {
      stop("Variable \"cuttrd\" must be of type \"numeric\"")
    }

    if (cuttrd > 1 | cuttrd < 0) {
      stop("Variable \"cuttrd\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(nagg) & !inherits(nagg, "numeric")) {
    stop("\"nagg\" must be of type \"numeric\"")
  }

  if (!is.null(startage) && !inherits(startage, "numeric")) {
    stop("\"startage\" must be of type \"numeric\"")
  }

  if (!is.null(ncase) & !inherits(ncase, "numeric")) {
    stop("\"ncase\" must be of type \"numeric\"")
  }
}
