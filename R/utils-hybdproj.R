#' Create cuttrend for hybdproj
#'
#' @inheritParams canproj
#' @param nonewpred Number of predicted periods
#'
#' @keywords internal
get_hybd_cuttrend <- function(shortp, nonewpred, nagg, cuttrd) {
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

#' Select hybrid model
#'
#' Select projection model for hybrid method
#'
#' @inheritParams canproj
#' @param casesB A `data.frame()`.
#' @param apdata A `data.frame()`.
#' @param dnoperiods Number of periods.
#'
#' @keywords internal
select_hybd_method <- function(casesB, apdata, nagg, pGOF, dnoperiods, pD) {
  if (ncol(casesB) == 5) {
    apdatan <- apdata
    lastper <- 5
    projbase <- 5 * nagg
  } else {
    mod1 <- stats::glm(
      Cases ~ as.factor(Age) + Period + offset(log(y)) - 1,
      data = apdata,
      family = stats::poisson
    )
    mod2 <- stats::glm(
      Cases ~ as.factor(Age) + as.factor(Age) * Period + offset(log(y)) - 1,
      data = apdata,
      family = stats::poisson
    )
    pdiff <- 1 -
      stats::pchisq(
        (mod1$deviance - mod2$deviance),
        (mod1$df.residual - mod2$df.residual)
      )
    if (is.null(pdiff)) {
      mod <- "common"
    } else {
      if (pdiff < pGOF) {
        mod <- "age-specific"
      } else {
        mod <- "common"
      }
    }

    trydat <- apdata
    likeh <- rep(0, (dnoperiods - 5))
    for (i in 1:(dnoperiods - 5)) {
      trydat$Period[trydat$Period <= i] <- 0
      if (mod == "common") {
        likeh[i] <- stats::glm(
          Cases ~ as.factor(Age) + Period + offset(log(y)) - 1,
          data = trydat,
          family = stats::poisson
        )$deviance /
          (-2)
      } else {
        likeh[i] <- stats::glm(
          Cases ~ as.factor(Age) + as.factor(Age) * Period + offset(log(y)) - 1,
          data = trydat,
          family = stats::poisson
        )$deviance /
          (-2)
      }
    }
    cuty <- which(likeh == max(likeh))
    apdatan <- apdata[apdata$Period >= cuty, ]
    if (cuty > 1) {
      apdatan$Period <- apdatan$Period - cuty + 1
    }
    lastper <- length(cuty:dnoperiods)
    projbase <- lastper * nagg
  }

  mode1 <- stats::glm(
    Cases ~ as.factor(Age) + Period + offset(log(y)) - 1,
    data = apdatan,
    family = stats::poisson
  )
  mode2 <- stats::glm(
    Cases ~ as.factor(Age) + as.factor(Age) * Period + offset(log(y)) - 1,
    data = apdatan,
    family = stats::poisson
  )
  pdiffn <- 1 -
    stats::pchisq(
      (mod1$deviance - mod2$deviance),
      (mod1$df.residual - mod2$df.residual)
    )
  pd1 <- summary(mode1)$coef["Period", 4]
  pd2 <- 1 - stats::pchisq(mode2$deviance, mode2$df.residual)

  if (is.null(pdiffn)) {
    if (pd1 > pD) {
      fmodel <- "average"
    } else {
      fmodel <- "common-trend"
    }
  } else {
    if (pdiffn > pGOF) {
      if (pd1 > pD) {
        fmodel <- "average"
      } else {
        fmodel <- "common-trend"
      }
    } else {
      if (pd2 > pD) {
        fmodel <- "age-specific"
      } else {
        fmodel <- "nba-specific"
      }
    }
  }

  return(list(fmodel, apdatan, lastper, cuty, projbase))
}
