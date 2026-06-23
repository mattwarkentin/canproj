#' Validate 'get prediction' inputs
#'
#' Check to see if inputs to 'get prediction' are valid
#'
#' @inheritParams canproj
#' @param nonewpred Number of predicted periods
#'
#' @keywords internal
get_hybd_cuttrend <- function(shortp, nonewpred, nagg) {
  cuttrend <- rep(shortp, nonewpred)

  if (nagg == 1) {
    if (nonewpred <= 5) {
      cuttrend <- rep(shortp, nonewpred)
    } else {
      for (i in 6:nonewpred) {
        cuttrend[i] <- shortp + (i - 5) * cuttrd
      }
    }
    cuttrend[cuttrend > 1] <- 1
  } else if (nagg == 2) {
    if (nonewpred <= 3) {
      cuttrend <- rep(shortp, nonewpred)
    } else {
      for (i in 4:nonewpred) {
        cuttrend[i] <- shortp + (i - 3) * 2 * cuttrd
      }
    }
    cuttrend[cuttrend > 1] <- 1
  } else if (nagg == 3) {
    if (nonewpred <= 2) {
      cuttrend <- rep(shortp, nonewpred)
    } else {
      for (i in 3:nonewpred) {
        cuttrend[i] <- shortp + (i - 2) * 3 * cuttrd
      }
    }
    cuttrend[cuttrend > 1] <- 1
  } else if (nagg == 4) {
    if (nonewpred == 1) {
      cuttrend <- shortp
    } else {
      for (i in 2:nonewpred) {
        cuttrend[i] <- shortp + (i - 1) * 4 * cuttrd
      }
    }
    cuttrend[cuttrend > 1] <- 1
  } else if (nagg == 5) {
    if (nonewpred == 1) {
      cuttrend <- shortp
    } else {
      for (i in 2:nonewpred) {
        cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
      }
    }
    cuttrend[cuttrend > 1] <- 1
  }

  return(cuttrend)
}
