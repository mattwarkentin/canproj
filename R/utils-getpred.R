#' Validate 'get prediction' inputs
#'
#' Check to see if inputs to 'get prediction' are valid
#'
#' @inheritParams canproj
#' @inheritParams hybdproj_get_predictions
#'
#' @keywords internal
validate_getpred_inputs <- function(byage, standpop, incidence, agegroups) {
  if ((!is.null(standpop)) && (!incidence)) {
    stop(
      "\"standpop\" should only be used with incidence predictions (incidence=T)"
    )
  }

  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)

    if (round(sum(standpop@weights), 5) != 1) {
      stop("\"standpop\" must be of sum 1")
    }

    if ((length(standpop) != length(agegroups)) && (agegroups[1] != "all")) {
      stop("\"standpop\" must be the same length as \"agegroups\"")
    }

    if (byage) {
      stop("\"standpop\" is only valid for \"byage=F\"")
    }
  }
}

#' Make prediction table
#'
#' Check to see if inputs to 'get prediction' are valid
#'
#' @inheritParams canproj
#' @inheritParams hybdproj_get_predictions
#' @param object "adpcproj", "hybdproj", or "acproj" object
#'
#' @keywords internal
make_pred_table <- function(
  object,
  agegroups,
  standpop,
  byage,
  incidence,
  excludeobs
) {
  datatable <- object$predictions
  pyr <- data.frame(object$pyr)

  if (agegroups[1] != "all") {
    datatable <- datatable[agegroups, ]
    pyr <- pyr[agegroups, ]
  }

  if (!is.null(standpop)) {
    incidence <- (datatable / pyr) * 100000
    incidence[is.na(incidence)] <- 0
    res <- colSums(incidence * standpop@weights)
  } else {
    if (!byage) {
      datatable <- colSums(datatable)
      pyr <- colSums(pyr)
    }

    if (incidence) {
      res <- (datatable / pyr) * 100000
      res[is.na(res)] <- 0
    } else {
      res <- datatable
    }
  }

  if (excludeobs) {
    predstart <- ncol(res) - object$nopred + 1
    res <- res[predstart:(predstart + object$nopred - 1)]
  }

  return(res)
}
