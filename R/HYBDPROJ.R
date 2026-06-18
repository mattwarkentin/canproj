#' HYBDPROJ
#'
#' R functions for projection of cancer incidence/mortality using the modified
#'   Hybrid methods. Modified Hybrid by adding choice of age-model, cut-trend
#'   parameter, and power 5 link function.
#'
#' @inheritParams canproj
#'
#' @return A `list()`.
#'
#' @export
hybdproj <- function(
  cdat,
  pdat,
  standpop,
  projfor = "incidence",
  nagg = NULL,
  ncase = NULL,
  cuttrd = 0.04,
  shortp = 0,
  linkfunc = "power5",
  pD = 0.05,
  pGOF = 0.05
) {
  S7::check_is_S7(standpop, StandardPopulation)

  if (is.null(ncase)) {
    if (projfor == "incidence") {
      ncase <- 1
    } else {
      ncase <- 3 / 5
    }
  }

  if (is.null(nagg)) {
    masr <- mean(obasr(cdat, pdat, standpop)[, 1])
    if (projfor == "incidence") {
      if (masr > 15) {
        nagg <- 1
      } else if (masr > 10) {
        nagg <- 2
      } else if (masr > 5) {
        nagg <- 3
      } else {
        nagg <- 4
      }
    } else {
      if (masr > 9) {
        nagg <- 1
      } else if (masr > 6) {
        nagg <- 2
      } else if (masr > 3) {
        nagg <- 3
      } else {
        nagg <- 4
      }
    }
  }

  aggdata <- datagg(cdat, pdat, nagg)
  cases <- aggdata$cases
  pyr <- aggdata$pyr

  percases <- ncol(cases)
  if (percases < 5) {
    stop("Minimum number of period is 5 (5*nagg years) in \"cases\"")
  }

  noperiod <- percases
  if (percases > 25) {
    noperiod <- 25
    cases <- cases[, -c(1:(percases - 25))]
    pyr <- pyr[, -c(1:(percases - 25))]
  }

  est <- hybdproj.estimate(
    cases = cases,
    pyr = pyr,
    nagg = nagg,
    ncase = ncase,
    linkfunc = linkfunc,
    pD = pD,
    pGOF = pGOF
  )
  pred <- hybdproj.prediction(
    hybdproj.estimate.object = est,
    cuttrd = cuttrd,
    shortp = shortp
  )
  return(pred)
}


#' hybdproj.estimate
#'
#' Define projection base, model selection, and fitting
#'
#' @inheritParams canproj
#' @param cases `data.frame` with number of cases in `nagg`-year period by ascending age groups in row.
#' @param pyr `data.frame` with observed and projected population size in `nagg`-year.
#'  period by ascending age groups in row.
#'
#' @return A `list()`.
#'
#' @export
hybdproj.estimate <- function(
  cases,
  pyr,
  nagg,
  ncase,
  linkfunc = "power5",
  pD = 0.05,
  pGOF = 0.05
) {
  if (dim(cases)[2] > dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")
  }

  if (dim(pyr)[2] == dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates")
  }

  if ((dim(pyr)[2] - dim(cases)[2]) > (30 / nagg)) {
    stop("Package can not project more than 30 years")
  }

  if (dim(cases)[2] < 5) {
    stop("\"noperiod\" must be 5 or larger")
  }

  nage <- 1:nrow(cases)
  mcase <- apply(cases, 1, mean)
  caseage <- as.data.frame(cbind(nage, cases))
  pyrage <- as.data.frame(cbind(nage, pyr))

  casesA <- caseage[mcase < (nagg * ncase), ]
  pyrA <- pyrage[mcase < (nagg * ncase), ]

  casesB <- caseage[mcase >= (nagg * ncase), ]
  pyrB <- pyrage[mcase >= (nagg * ncase), ]

  dnoperiods <- ncol(casesB) - 1
  dnoagegr <- nrow(casesB)
  agpreg <- casesB[, 1]
  agpave <- casesA[, 1]

  if (dnoagegr < 2) {
    stop("\"dnoagegr\" must be 2 or larger")
  }

  ageno <- rep(agpreg, dnoperiods)
  periodno <- sort(rep(1:dnoperiods, dnoagegr))
  y <- c(as.matrix(pyrB[, 2:(dnoperiods + 1)]))
  apdata <- data.frame(
    Age = ageno,
    Period = periodno,
    Cases = c(as.matrix(casesB[, -1])),
    y = y
  )

  if (dim(casesB)[2] == 5) {
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

  if (fmodel == "common-trend") {
    if (linkfunc == "power5") {
      y <- apdatan$y
      power5link <- stats::poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu) {
        (mu / y)^0.2
      }
      power5link$linkinv <- function(eta) {
        pmax(.Machine$double.eps, y * eta^5)
      }
      power5link$mu.eta <- function(eta) {
        pmax(.Machine$double.eps, 5 * y * eta^4)
      }

      res.glm <- stats::glm(
        Cases ~ as.factor(Age) + Period - 1,
        data = apdatan,
        family = power5link
      )
    } else if (linkfunc == "log") {
      res.glm <- stats::glm(
        Cases ~ as.factor(Age) + Period + offset(log(y)) - 1,
        data = apdatan,
        family = stats::poisson(link = log)
      )
    } else if (linkfunc == "sqrt") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) + Period - 1,
          data = apdatan,
          family = stats::poisson(link = sqrt)
        )
      )
    } else if (linkfunc == "identity") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) + Period - 1,
          data = apdatan,
          family = stats::poisson(link = identity)
        )
      )
    } else {
      stop("Unknown \"linkfunc\"")
    }

    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  } else if (fmodel == "age-specific") {
    if (linkfunc == "power5") {
      y <- apdatan$y
      power5link <- stats::poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu) {
        (mu / y)^0.2
      }
      power5link$linkinv <- function(eta) {
        pmax(.Machine$double.eps, y * eta^5)
      }
      power5link$mu.eta <- function(eta) {
        pmax(.Machine$double.eps, 5 * y * eta^4)
      }

      res.glm <- stats::glm(
        Cases ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
        data = apdatan,
        family = power5link
      )
    } else if (linkfunc == "log") {
      res.glm <- stats::glm(
        Cases ~
          as.factor(Age) +
            as.factor(Age) * Period +
            offset(log(y)) -
            1 -
            Period,
        data = apdatan,
        family = stats::poisson(link = log)
      )
    } else if (linkfunc == "sqrt") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
          data = apdatan,
          family = stats::poisson(link = sqrt)
        )
      )
    } else if (linkfunc == "identity") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
          data = apdatan,
          family = stats::poisson(link = identity)
        )
      )
    } else {
      stop("Unknown \"linkfunc\"")
    }

    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  } else if (fmodel == "nba-specific") {
    if (linkfunc == "power5") {
      y <- apdatan$y
      glmnb <- suppressWarnings(MASS::glm.nb(
        Cases ~
          as.factor(Age) +
            as.factor(Age) * Period -
            1 -
            Period +
            offset(log(y)),
        data = apdatan,
        link = log
      ))
      theta <- as.numeric(MASS::theta.md(
        apdatan$Cases,
        stats::fitted(glmnb),
        dfr = stats::df.residual(glmnb)
      ))
      nbpower5link <- MASS::negative.binomial(theta)
      nbpower5link$link <- "0.2 root link negative.binomial(theta) family"
      nbpower5link$linkfun <- function(mu) {
        (mu / y)^0.2
      }
      nbpower5link$linkinv <- function(eta) {
        pmax(.Machine$double.eps, y * eta^5)
      }
      nbpower5link$mu.eta <- function(eta) {
        pmax(.Machine$double.eps, 5 * y * eta^4)
      }

      res.glm <- stats::glm(
        Cases ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
        data = apdatan,
        family = nbpower5link
      )
    } else if (linkfunc == "log") {
      res.glm <- suppressWarnings(MASS::glm.nb(
        Cases ~
          as.factor(Age) +
            as.factor(Age) * Period +
            offset(log(y)) -
            1 -
            Period,
        data = apdatan,
        link = log
      ))
    } else if (linkfunc == "sqrt") {
      res.glm <- suppressWarnings(MASS::glm.nb(
        Cases / y ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
        data = apdatan,
        link = sqrt
      ))
    } else if (linkfunc == "identity") {
      res.glm <- suppressWarnings(MASS::glm.nb(
        Cases / y ~ as.factor(Age) + as.factor(Age) * Period - 1 - Period,
        data = apdatan,
        link = identity
      ))
    } else {
      stop("Unknown \"linkfunc\"")
    }
    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  } else {
    if (linkfunc == "power5") {
      y <- apdatan$y
      power5link <- stats::poisson()
      power5link$link <- "0.2 root link Poisson family"
      power5link$linkfun <- function(mu) {
        (mu / y)^0.2
      }
      power5link$linkinv <- function(eta) {
        pmax(.Machine$double.eps, y * eta^5)
      }
      power5link$mu.eta <- function(eta) {
        pmax(.Machine$double.eps, 5 * y * eta^4)
      }
      #
      res.glm <- stats::glm(
        Cases ~ as.factor(Age) - 1,
        data = apdatan,
        family = power5link
      )
    } else if (linkfunc == "log") {
      res.glm <- stats::glm(
        Cases ~ as.factor(Age) + offset(log(y)) - 1,
        data = apdatan,
        family = stats::poisson(link = log)
      )
    } else if (linkfunc == "sqrt") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) - 1,
          data = apdatan,
          family = stats::poisson(link = sqrt)
        )
      )
    } else if (linkfunc == "identity") {
      suppressWarnings(
        res.glm <- stats::glm(
          Cases / y ~ as.factor(Age) - 1,
          data = apdatan,
          family = stats::poisson(link = identity)
        )
      )
    } else {
      stop("Unknown \"linkfunc\"")
    }
    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  }

  res <- list(
    glm = res.glm,
    cases = cases,
    pyr = pyr,
    agrpave = agpave,
    lastper = lastper,
    cuty = cuty,
    noperiod = lastper,
    noyearagg = nagg,
    nocaseagp = ncase,
    linkfunc = linkfunc,
    agrpmod = agpreg,
    projbase = projbase,
    finalmod = fmodel,
    gofpvalue = pvalue
  )
  class(res) <- "hybdproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' hybdproj.prediction
#'
#' Extrapolate estimated trend from the final model
#'
#' @inheritParams canproj
#'
#' @param hybdproj.estimate.object An object based on the `hybdproj.estimate` function.
#'
#' @return A `list()`.
#'
#' @export
hybdproj.prediction <- function(
  hybdproj.estimate.object,
  cuttrd = 0.04,
  shortp = 0
) {
  if (!inherits(hybdproj.estimate.object, "hybdproj.estimate")) {
    stop(
      "Variable \"hybdproj.estimate.object\" must be of type \"hybdproj.estimate\""
    )
  }

  cases <- hybdproj.estimate.object$cases
  pyr <- hybdproj.estimate.object$pyr
  nagg <- hybdproj.estimate.object$noyearagg
  ncase <- hybdproj.estimate.object$nocaseagp
  projbase <- hybdproj.estimate.object$projbase
  agpreg <- hybdproj.estimate.object$agrpmod
  agpave <- hybdproj.estimate.object$agrpave
  fmodel <- hybdproj.estimate.object$finalmod
  noperiod <- hybdproj.estimate.object$noperiod
  lastper <- hybdproj.estimate.object$lastper
  ngroups <- nrow(cases)

  noobsper <- ncol(cases)
  nototper <- ncol(pyr)
  nonewpred <- nototper - noobsper

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

  driftmp <- cumsum(1 - cuttrend)

  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode", 1:nototper)
    } else {
      years <- dimnames(pyr)[[2]]
    }
  }

  datatable <- data.frame(cases, matrix(NA, nrow(cases), nonewpred))
  names(datatable) <- years

  obsrat <- data.frame(matrix(0, dim(cases)[1], 5))

  for (age in 1:(dim(cases)[1])) {
    obsrat[age, ] <- as.matrix(cases[age, (noobsper - 4):noobsper]) /
      as.matrix(pyr[age, (noobsper - 4):noobsper])
  }
  obsrate <- apply(obsrat, 1, mean)

  if (fmodel == "age-specific" || fmodel == "nba-specific") {
    coef <- hybdproj.estimate.object$glm$coef
    m.eff <- cbind(
      coef[1:(length(coef) / 2)],
      coef[(length(coef) / 2 + 1):length(coef)],
      agpreg
    )
    row.names(m.eff) <- NULL
    avef <- cbind(rep(NA, length(agpave)), rep(0, length(agpave)), agpave)
    mcoef <- data.frame(rbind(avef, m.eff))
    acoef <- mcoef[with(mcoef, order(mcoef[, 3])), ]
    row.names(acoef) <- NULL
    colnames(acoef) <- c("a.eff", "p.eff", "agrp")

    for (age in 1:19) {
      if (is.na(acoef$a.eff[acoef$agrp == age])) {
        rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc == "power5") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^5
        } else if (hybdproj.estimate.object$linkfunc == "log") {
          rate <- exp(acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        } else if (hybdproj.estimate.object$linkfunc == "sqrt") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^2
        } else {
          # identity link:
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        }
      }
      datatable[age, (noobsper + 1):nototper] <- rate *
        pyr[age, (noobsper + 1):nototper]
    }
  } else if (fmodel == "common-trend") {
    coef <- hybdproj.estimate.object$glm$coef
    m.eff <- cbind(
      coef[1:(length(coef) - 1)],
      rep(coef[length(coef)], length(agpreg)),
      agpreg
    )
    row.names(m.eff) <- NULL
    avef <- cbind(rep(NA, length(agpave)), rep(0, length(agpave)), agpave)
    mcoef <- data.frame(rbind(avef, m.eff))
    acoef <- mcoef[with(mcoef, order(mcoef[, 3])), ]
    row.names(acoef) <- NULL
    colnames(acoef) <- c("a.eff", "p.eff", "agrp")

    for (age in 1:nrow(cases)) {
      if (is.na(acoef$a.eff[acoef$agrp == age])) {
        rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc == "power5") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^5
        } else if (hybdproj.estimate.object$linkfunc == "log") {
          rate <- exp(acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        } else if (hybdproj.estimate.object$linkfunc == "sqrt") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^2
        } else {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        }
      }
      datatable[age, (noobsper + 1):nototper] <- rate *
        pyr[age, (noobsper + 1):nototper]
    }
  } else {
    coef <- hybdproj.estimate.object$glm$coef
    m.eff <- cbind(coef[1:length(coef)], rep(0, length(agpreg)), agpreg)
    row.names(m.eff) <- NULL
    avef <- cbind(rep(NA, length(agpave)), rep(0, length(agpave)), agpave)
    mcoef <- data.frame(rbind(avef, m.eff))
    acoef <- mcoef[with(mcoef, order(mcoef[, 3])), ]
    row.names(acoef) <- NULL
    colnames(acoef) <- c("a.eff", "p.eff", "agrp")

    for (age in 1:nrow(cases)) {
      if (is.na(acoef$a.eff[acoef$agrp == age])) {
        rate <- rep(obsrate[age], length(driftmp))
      } else {
        if (hybdproj.estimate.object$linkfunc == "power5") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^5
        } else if (hybdproj.estimate.object$linkfunc == "log") {
          rate <- exp(acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        } else if (hybdproj.estimate.object$linkfunc == "sqrt") {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))^2
        } else {
          rate <- (acoef$a.eff[age] + acoef$p.eff[age] * (lastper + driftmp))
        }
      }
      datatable[age, (noobsper + 1):nototper] <- rate *
        pyr[age, (noobsper + 1):nototper]
    }
  }

  res <- list(
    predictions = datatable,
    pyr = pyr,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    nopred = nonewpred,
    noperiod = hybdproj.estimate.object$noperiod,
    lastperiod = lastper,
    noobsper = noobsper,
    nototper = nototper,
    noyearagg = hybdproj.estimate.object$noyearagg,
    nocaseagp = hybdproj.estimate.object$nocaseagp,
    agrpave = hybdproj.estimate.object$agrpave,
    agrpmod = hybdproj.estimate.object$agrpmod,
    linkfunc = hybdproj.estimate.object$linkfunc,
    projbase = hybdproj.estimate.object$projbase,
    finalmod = hybdproj.estimate.object$finalmod,
    gofpvalue = hybdproj.estimate.object$gofpvalue,
    glm = hybdproj.estimate.object$glm
  )

  class(res) <- "hybdproj"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' hybdproj.getpred
#'
#' Extract projection results by n-year period
#'
#' @inheritParams canproj
#' @param incidence Whether to give rates (`T`) or numbers (`F`).
#' @param hybdproj.object An object based on the 'hybdproj()' function.
#' @param excludeobs Whether to include observed values (`T` or `F`).
#' @param byage Report numbers by age groups (`T`), or use age-standardized rates (`F`).
#' @param agegroups Age groups to include. "all" (default) includes all age groups,
#'  and individual groups can be selected by group number. E.g. c(1:3, 7) includes
#'  the first 3 groups and the seventh group.
#'
#' @return A `data.frame()`.
#'
#' @export
hybdproj.getpred <- function(
  hybdproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
) {
  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  if (!inherits(hybdproj.object, "hybdproj")) {
    stop("Variable \"hybdproj.object\" must be of type \"hybdproj\"")
  }

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

  datatable <- hybdproj.object$predictions
  pyr <- data.frame(hybdproj.object$pyr)

  if (agegroups[1] != "all") {
    datatable <- datatable[agegroups, ]
    pyr <- pyr[agegroups, ]
  }

  if (!is.null(standpop)) {
    datainc <- (datatable / pyr) * 100000
    datainc[is.na(datainc)] <- 0
    res <- apply(datainc * standpop@weights, 2, sum)
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
    if (is.matrix(res)) {
      predstart <- ncol(res) - hybdproj.object$nopred + 1
      res <- res[, predstart:(predstart + hybdproj.object$nopred - 1)]
    } else {
      predstart <- length(res) - hybdproj.object$nopred + 1
      res <- res[predstart:(predstart + hybdproj.object$nopred - 1)]
    }
  }

  return(res)
}


#' hybdproj.getproj
#'
#' Extract projection results.
#'
#' @inheritParams hybdproj.getpred
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @export
hybdproj.getproj <- function(
  cdat,
  pdat,
  startp,
  hybdproj.object,
  standpop = NULL,
  Ave5 = NULL,
  sum5 = NULL
) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  finalmod <- hybdproj.object$finalmod

  r0 <- hybdproj.getpred(hybdproj.object, incidence = T)

  nagg <- hybdproj.object$noyearagg

  outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = nagg)

  if (finalmod != "average") {
    if (is.null(standpop)) {
      return(outasp)
    } else {
      outann <- asry(outasp, pdat, standpop)
      return(outann)
    }
  } else {
    if (is.null(Ave5)) {
      if (is.null(standpop)) {
        return(outasp)
      } else {
        outann <- asry(outasp, pdat, standpop)
        return(outann)
      }
    } else {
      mod <- ave5proj(cdat, pdat, startp, sum5 = sum5)
      outasp <- mod$agsproj
      if (is.null(standpop)) {
        return(outasp)
      } else {
        outann <- asry(outasp, pdat, standpop)
        return(outann)
      }
    }
  }
}


#' summary.hybdproj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the 'hybdproj()' function.
#' @param printpred Whether to print the observed and predicted number of cases (`T` or `F`).
#' @param printcall Whether to print function `Call` for hybdproj.object (`T` or `F`).
#' @param digits Number of digits in output, default (`0`) is integer only.
#' @param ... Other parameters.
#'
#' @return An information table describing the method used.
#'
#' @export
summary.hybdproj <- function(
  object,
  printpred = F,
  printcall = F,
  digits = 0,
  ...
) {
  method <- "Hybrid approach"

  if (!inherits(object, "hybdproj")) {
    stop("Variable \"object\" must be of type \"hybdproj\"")
  }

  nototper <- object$nototper
  noobsper <- object$noobsper
  nopred <- object$nopred
  noypred <- (object$nopred) * (object$noyearagg)
  noperiod <- object$noperiod
  obsto <- names(object$predictions)[
    dim(object$predictions)[2] - nopred
  ]
  predcases <- object$predictions[, (noobsper + 1):(nototper)]

  if (!is.null(object$gofpvalue)) {
    gofpvalue <- round(object$gofpvalue, 4)
  } else {
    gofpvalue <- NA
  }

  if (printpred) {
    cat("Predicted number of cases or deaths:")
    cat("(observations up to", obsto, ")\n")
    print(round(as.matrix(predcases), digits = digits))
    cat("\n")
  }
  cat("\nPrediction done with:\n")

  moptions <- matrix(NA, 12, 2)

  moptions[, 1] <- c(
    "Method:",
    "Number of prediction years:",
    "First period cutting trend:",
    "Degenerating trend per year:",
    "Projection base (years):",
    "Aggregating years (nagg):",
    "Age-cases per year (ncase):",
    "Model for regression:",
    "Link function for GLM:",
    "P-value for goodness of fit:",
    "Age group for regression:",
    "Age group for average method:"
  )

  moptions[, 2] <- c(
    method,
    noypred,
    object$shortp,
    object$cuttrd,
    object$projbase,
    object$noyearagg,
    object$nocaseagp,
    object$finalmod,
    object$linkfunc,
    gofpvalue,
    paste(object$agrpmod, collapse = ","),
    paste(object$agrpave, collapse = ",")
  )

  maxl <- max(nchar(moptions[, 1]))

  for (i in 1:dim(moptions)[1]) {
    spaces <- rep(" ", maxl - nchar(moptions[i, 1]) + 2)
    cat(moptions[i, 1], spaces, moptions[i, 2], "\n", sep = "")
  }

  if (printcall) {
    cat("\n  Call: ")
    dput(attr(object, "Call"))
  }

  invisible(object)
}


#' glm.hybdproj
#'
#' Summarize estimations from the final model.
#'
#' @inheritParams hybdproj.getpred
#'
#' @return A summary table from the `glm.object`.
#'
#' @export
glm.hybdproj <- function(hybdproj.object) {
  summary(hybdproj.object$glm)
}


#' plot.hybdproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' hybdproj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
#' @param x An object based on the 'hybdproj()' function.
#'
#' @export
plot.hybdproj <- function(
  x,
  cdat,
  pdat,
  startp,
  standpop,
  startplot = 1,
  xlab = "Calendar Year",
  ylab = "Rates",
  main = "",
  labels = NULL,
  ylim = NULL,
  lty = c(1, 3),
  col = c(1, 1),
  new = T,
  ...
) {
  if (!inherits(x, "hybdproj")) {
    stop("Variable \"x\" must be of type \"hybdproj\"")
  }
  S7::check_is_S7(standpop, StandardPopulation)

  indat <- hybdproj.getproj(
    cdat,
    pdat,
    startp = startp,
    x,
    standpop = standpop
  )
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  obsy <- dim(cdat)[2]
  nopredy <- length(indata) - obsy
  if (is.null(labels)) {
    labels <- row.names(indat)
  }

  maxx <- length(indata)
  if (new) {
    maxy <- max(indata) * (21 / 20)
    if (is.null(ylim)) {
      ylim <- c(0, maxy)
    }
    graphics::plot(
      c(1, maxx),
      ylim,
      type = "n",
      ylab = ylab,
      xlab = xlab,
      axes = F,
      ...
    )
    graphics::axis(2)
    graphics::axis(1, at = 1:maxx, labels = labels)
    graphics::box()
    graphics::title(main)
  }

  graphics::lines(
    1:(maxx - nopredy),
    indata[1:(maxx - nopredy)],
    type = "o",
    pch = 20,
    lty = lty[1],
    col = col[1],
    ...
  )

  graphics::lines(
    (maxx - nopredy):maxx,
    indata[(maxx - nopredy):maxx],
    lty = lty[2],
    col = col[2],
    ...
  )

  invisible(x)
}
