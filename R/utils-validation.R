#' Validate inputs for summary functions
#'
#' @keywords internal
validate_summary_inputs <- function(
  printpred = NULL,
  printcall = NULL,
  digits = NULL
) {
  if (!is.null(printpred) && !inherits(printpred, "logical")) {
    rlang::abort("Variable \"printpred\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!is.null(printcall) && !inherits(printcall, "logical")) {
    rlang::abort("Variable \"printcall\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!is.null(digits)) {
    if (!inherits(digits, "numeric")) {
      rlang::abort("Variable \"digits\" must be of type \"numeric\"")
    }

    if (digits < 0) {
      rlang::abort("Variable \"digits\" must be >= 0")
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
      rlang::abort(
        "Variable \"standpop\" must have same number of age groups as \"pdat\""
      )
    }
  }

  if (!inherits(pdat, "data.frame") & !inherits(pdat, "matrix")) {
    rlang::abort("\"pdat\" must be of type \"data.frame\" or \"matrix\"")
  }

  if (!is.null(object$predictions)) {
    if (nrow(pdat) != nrow(object$predictions)) {
      rlang::abort(
        "\"pdat\" must have the same number of age groups as predictions"
      )
    }
  } else {
    if (nrow(pdat) != nrow(object$agsproj)) {
      rlang::abort(
        "\"pdat\" must have same number of age groups as predictions"
      )
    }
  }

  if (!is.null(cdat)) {
    if (!inherits(cdat, "data.frame") & !inherits(cdat, "matrix")) {
      rlang::abort("\"cdat\" must be of type \"data.frame\" or \"matrix\"")
    }

    if (nrow(cdat) != nrow(object$predictions)) {
      rlang::abort(
        "\"cdat\" must have the same number of age groups as predictions"
      )
    }
  }

  if (!is.null(startp)) {
    if (!inherits(startp, "numeric")) {
      rlang::abort("\"startp\" must be of type \"numeric\"")
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
    rlang::abort("Variable \"incidence\" must be \"TRUE\" or \"FALSE\"")
  }

  if (!inherits(excludeobs, "logical")) {
    rlang::abort("Variable \"excludeobs\" must be \"TRUE\" or \"FALSE\"")
  }

  if ((!is.null(standpop)) && (!incidence)) {
    rlang::abort(
      "\"standpop\" should only be used with incidence predictions (incidence=T)"
    )
  }

  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)

    if ((length(standpop) != length(agegroups)) && (agegroups[1] != "all")) {
      rlang::abort("\"standpop\" must be the same length as \"agegroups\"")
    }

    if (byage) {
      rlang::abort("\"standpop\" is only valid for \"byage=F\"")
    }
  }
}

#' Validate inputs for prediction functions
#'
#' @keywords internal
validate_prediction_inputs <- function(cuttrd, shortp) {
  if (!inherits(cuttrd, "numeric")) {
    rlang::abort("Variable \"cuttrd\" must be of type \"numeric\"")
  }

  if (cuttrd > 1 | cuttrd < 0) {
    rlang::abort("Variable \"cuttrd\" must be >= 0 and <= 1")
  }

  if (!inherits(shortp, "numeric")) {
    rlang::abort("Variable \"shortp\" must be of type \"numeric\"")
  }

  if (shortp > 1 | shortp < 0) {
    rlang::abort("Variable \"shortp\" must be >= 0 and <= 1")
  }
}

#' Validate inputs for estimate functions
#'
#' @keywords internal
validate_estimate_inputs <- function(pyr, cases, pGOF, linkfunc) {
  if (ncol(cases) > ncol(pyr)) {
    rlang::abort(
      "\"pyr\" must include information about all periods in \"cases\""
    )
  }

  if (ncol(pyr) == ncol(cases)) {
    rlang::abort("\"pyr\" must include information on future rates")
  }

  if (nrow(cases) != nrow(pyr)) {
    rlang::abort(
      "\"cases\" and \"pyr\" must have the same number of age groups"
    )
  }

  if (ncol(cases) < 2) {
    rlang::abort("\"cases\" must have at least 2 age groups")
  }

  if (ncol(pyr) < 2) {
    rlang::abort("\"pyr\" must have at least 2 age groups")
  }

  if (!inherits(pGOF, "numeric")) {
    rlang::abort("Variable \"pGOF\" must be of type \"numeric\"")
  }

  if (pGOF > 1 | pGOF < 0) {
    rlang::abort("Variable \"pGOF\" must be >= 0 and <= 1")
  }

  if (!(linkfunc) %in% list("power5", "sqrt", "log", "identity")) {
    rlang::abort(
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

  if (!is.null(n5case) && !inherits(n5case, "numeric")) {
    rlang::abort("\"n5case\" must be of type \"numeric\" or NULL")
  }

  if (!is.null(projfor) && !(projfor %in% list("incidence", "mortality"))) {
    rlang::abort("\"projfor\" must be one of \"incidence\" or \"mortality\"")
  }

  if (
    !is.null(linkfunc) &&
      !(linkfunc) %in% list("power5", "sqrt", "log", "identity")
  ) {
    rlang::abort(
      "Variable \"linkfunc\" must be one of \"power5\", \"log\", \"sqrt\", or \"identity\""
    )
  }

  if (!is.null(startp) && !inherits(startp, "numeric")) {
    rlang::abort("\"startp\" must be of type \"numeric\"")
  }

  if (!is.null(newcohort) && !inherits(newcohort, "logical")) {
    rlang::abort("\"newcohort\" must be of type \"logical\"")
  }

  if (!is.null(pGOF)) {
    if (!inherits(pGOF, "numeric")) {
      rlang::abort("Variable \"pGOF\" must be of type \"numeric\"")
    }

    if (pGOF > 1 | pGOF < 0) {
      rlang::abort("Variable \"pGOF\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(pD)) {
    if (!inherits(pD, "numeric")) {
      rlang::abort("Variable \"pD\" must be of type \"numeric\"")
    }

    if (pD > 1 | pD < 0) {
      rlang::abort("Variable \"pD\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(shortp)) {
    if (!inherits(shortp, "numeric")) {
      rlang::abort("Variable \"shortp\" must be of type \"numeric\"")
    }

    if (shortp > 1 | shortp < 0) {
      rlang::abort("Variable \"shortp\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(cuttrd)) {
    if (!inherits(cuttrd, "numeric")) {
      rlang::abort("Variable \"cuttrd\" must be of type \"numeric\"")
    }

    if (cuttrd > 1 | cuttrd < 0) {
      rlang::abort("Variable \"cuttrd\" must be >= 0 and <= 1")
    }
  }

  if (!is.null(nagg) & !inherits(nagg, "numeric")) {
    rlang::abort("\"nagg\" must be of type \"numeric\"")
  }

  if (!is.null(startage) && !inherits(startage, "numeric")) {
    rlang::abort("\"startage\" must be of type \"numeric\"")
  }

  if (!is.null(ncase) & !inherits(ncase, "numeric")) {
    rlang::abort("\"ncase\" must be of type \"numeric\"")
  }
}
