#' ACPROJ
#'
#' R functions for projection of cancer incidence/mortality. Revising and
#'   combining nordpred and Osmond's to extrapolation cohort when no drift
#'   appears from nordpred.
#'
#' @inheritParams canproj
#' @param n5case Minimum number of cancer cases/deaths per 5 years for splitting data.
#'
#' @return A `list()`.
#'
#' @export
acproj <- function(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  startestage = NULL,
  cuttrd = 0.04,
  shortp = 0,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  if (is.null(n5case)) {
    if (projfor == "incidence") {
      n5case <- 5
    } else {
      n5case <- 3
    }
  }

  aggdata <- datagg(cdat, pdat, 5)
  cases <- aggdata$cases
  pyr <- aggdata$pyr

  if (is.null(startestage)) {
    dat <- as.matrix(cases)
    iage <- 1
    while (mean(dat[iage, ]) < n5case | dat[iage, 1] == 0) {
      iage <- iage + 1
    }
    startestage <- iage
  }

  percases <- ncol(cases)
  if (percases < 3) {
    stop("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  est <- acproj.estimate(
    cases = cases,
    pyr = pyr,
    noperiod = percases,
    pGOF = pGOF,
    startestage = startestage,
    linkfunc = linkfunc
  )

  pred <- acproj.prediction(
    acproj.estimate.object = est,
    cuttrd = cuttrd,
    shortp = shortp
  )

  return(pred)
}


#' acproj.estimate
#'
#' Fit age-cohort models
#'
#' @inheritParams canproj
#' @param cases `data.frame` with number of cases in `nagg`-year period by ascending age groups in row.
#' @param pyr `data.frame` with observed and projected population size in `nagg`-year.
#' @param noperiod Number of 5-year periods in historical data.
#'
#' @return A `list()`.
#'
#' @export
acproj.estimate <- function(
  cases,
  pyr,
  noperiod,
  startestage,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  if (dim(cases)[2] > dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")
  }

  if (dim(pyr)[2] == dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates")
  }

  if ((dim(pyr)[2] - dim(cases)[2]) > 6) {
    stop("Package can not project more than 6 periods")
  }

  if ((dim(cases)[2] - noperiod) < 0) {
    stop("More periods specified in \"noperiod\" than columns in \"cases\"")
  }

  if (noperiod < 3) {
    stop("\"noperiod\" must be 3 or larger")
  }

  dnoperiods <- ncol(cases)
  dnoagegr <- nrow(cases)

  ageno <- rep(1:dnoagegr, dnoperiods)
  periodno <- sort(rep(1:dnoperiods, dnoagegr))
  cohort <- max(ageno) - ageno + periodno
  y <- c(as.matrix(pyr[, 1:dnoperiods]))
  apcdata <- data.frame(
    Age = ageno,
    Cohort = cohort,
    Period = periodno,
    Cases = c(as.matrix(cases)),
    y = y
  )

  ## Selecting data for regression:
  apcdata <- apcdata[apcdata$Age >= startestage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]
  maxc <- max(apcdata$Cohort)
  midc <- ceiling(maxc / 2)

  if (linkfunc == "power5") {
    y <- apcdata$y

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
      Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
      data = apcdata,
      family = power5link
    )
  } else if (linkfunc == "log") {
    res.glm <- stats::glm(
      Cases ~ factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
      data = apcdata,
      family = stats::poisson(link = log)
    )
  } else if (linkfunc == "sqrt") {
    suppressWarnings(
      res.glm <- stats::glm(
        Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        family = stats::poisson(link = sqrt)
      )
    )
  } else if (linkfunc == "identity") {
    suppressWarnings(
      res.glm <- stats::glm(
        Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        family = stats::poisson(link = identity)
      )
    )
  } else {
    stop("Unknown \"linkfunc\"")
  }

  pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

  if (pvalue < pGOF) {
    distn <- "Negative-binomial"
    if (linkfunc == "power5") {
      y <- apcdata$y

      suppressWarnings(
        glmnb <- MASS::glm.nb(
          Cases ~
            factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
          data = apcdata,
          link = log
        )
      )

      theta <- as.numeric(MASS::theta.md(
        apcdata$Cases,
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
        Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        family = nbpower5link
      )
    } else if (linkfunc == "log") {
      suppressWarnings(
        res.glm <- MASS::glm.nb(
          Cases ~
            factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
          data = apcdata,
          link = log
        )
      )
    } else if (linkfunc == "sqrt") {
      suppressWarnings(
        res.glm <- MASS::glm.nb(
          Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
          data = apcdata,
          link = sqrt
        )
      )
    } else if (linkfunc == "identity") {
      suppressWarnings(
        res.glm <- MASS::glm.nb(
          Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
          data = apcdata,
          link = identity
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
    maxc = maxc,
    midc = midc,
    noperiod = noperiod,
    linkfunc = linkfunc,
    startestage = startestage,
    distribution = distn,
    gofpvalue = pvalue
  )
  class(res) <- "acproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' acproj.prediction
#'
#' Extrapolate estimated trend from age-cohort model
#'
#' @inheritParams canproj
#' @param acproj.estimate.object An object based on the `acproj.estimate()` function.
#'
#' @return A `list()`.
#'
#' @export
acproj.prediction <- function(
  acproj.estimate.object,
  cuttrd = 0.04,
  shortp = 0
) {
  if (!inherits(acproj.estimate.object, "acproj.estimate")) {
    stop(
      "Variable \"acproj.estimate.object\" must be of type \"acproj.estimate\""
    )
  }

  cases <- acproj.estimate.object$cases
  pyr <- acproj.estimate.object$pyr
  noperiod <- acproj.estimate.object$noperiod
  startestage <- acproj.estimate.object$startestage
  maxc <- acproj.estimate.object$maxc
  midc <- acproj.estimate.object$midc
  nototper <- ncol(pyr)
  noobsper <- ncol(cases)
  nonewpred <- nototper - noobsper
  ngroups <- nrow(cases)

  cuttrend <- rep(shortp, nonewpred)
  for (i in 1:nonewpred) {
    cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
  }
  cuttrend[cuttrend > 1] <- 1

  datatable <- data.frame(cases, matrix(NA, ngroups, nonewpred))
  names(datatable) <- names(pyr)

  skipped_ages <- 1:(startestage - 1)
  obsinc <- cases[skipped_ages, (noobsper - 1):noobsper] /
    pyr[skipped_ages, (noobsper - 1):noobsper]
  obsinc[is.na(obsinc)] <- 0

  datatable[skipped_ages, (noobsper + 1):nototper] <- ((obsinc[, 1] +
    obsinc[, 2]) /
    2) *
    pyr[skipped_ages, (noobsper + 1):nototper]

  coefficients <- acproj.estimate.object$glm$coefficients

  num_ages <- ngroups - startestage + 1
  ageff <- coefficients[1:num_ages]

  coeff <- c(
    coefficients[(num_ages + 1):(num_ages + midc - 1)],
    0,
    coefficients[(num_ages + midc):length(coefficients)]
  )
  lncoh <- length(coeff)

  # For new young cohorts in "nonewpred" by linear regression without the first and last 4 cohorts:
  ceff <- coeff[-c(1, 2, lncoh - 1, lncoh)]
  lnceff <- length(ceff)
  tt <- 1:lnceff
  coeflm <- stats::coef(stats::lm(ceff ~ tt - 1))
  ncoeff <- coeflm *
    c((lncoh - 3):(lncoh - 2), (cumsum(1 - cuttrend) + lncoh - 2))
  coheff <- c(coeff[1:(lncoh - 2)], ncoeff)

  for (age in startestage:ngroups) {
    agepar <- as.numeric(ageff[age - startestage + 1])
    coh <- (ngroups - startestage) -
      (age - startestage) +
      (noperiod + 1:nonewpred)
    cohpar <- coheff[coh]

    if (acproj.estimate.object$linkfunc == "power5") {
      rate <- (agepar + cohpar)^5
    } else if (acproj.estimate.object$linkfunc == "log") {
      rate <- exp(agepar + cohpar)
    } else if (acproj.estimate.object$linkfunc == "sqrt") {
      rate <- (agepar + cohpar)^2
    } else {
      rate <- (agepar + cohpar)
    }
    datatable[age, (noobsper + 1):nototper] <- rate *
      pyr[age, (noobsper + 1):nototper]
  }

  res <- list(
    predictions = datatable,
    pyr = pyr,
    nopred = nonewpred,
    noperiod = acproj.estimate.object$noperiod,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    gofpvalue = acproj.estimate.object$gofpvalue,
    distribution = acproj.estimate.object$distribution,
    startestage = acproj.estimate.object$startestage,
    glm = acproj.estimate.object$glm
  )
  class(res) <- "acproj"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' acproj.getpred
#'
#' Extract projection results by 5-year period
#'
#' @inheritParams canproj
#' @param acproj.object An object based on the `acproj()` function.
#' @param incidence Whether to give rates (`T`) or numbers (`F`).
#' @param excludeobs Whether to include observed values (`T` or `F`).
#' @param byage Report numbers by age groups (`T`), or use age-standardized rates (`F`).
#' @param agegroups Age groups to include. "all" (default) includes all age groups,
#'  and individual groups can be selected by group number. E.g. c(1:3, 7) includes
#'  the first 3 groups and the seventh group.
#'
#' @return A `data.frame()`.
#'
#' @export
acproj.getpred <- function(
  acproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
) {
  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  if (!inherits(acproj.object, "acproj")) {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")
  }

  validate_getpred_inputs(byage, standpop, incidence, agegroups)

  res <- make_pred_table(
    acproj.object,
    agegroups,
    standpop,
    byage,
    incidence,
    excludeobs
  )

  return(res)
}


#' acproj.getproj
#'
#' Extract annual projection results
#'
#' @inheritParams canproj
#' @inheritParams acproj.getpred
#'
#' @return A `data.frame()`.
#'
#' @export
acproj.getproj <- function(cdat, pdat, startp, acproj.object, standpop = NULL) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  r0 <- acproj.getpred(acproj.object, incidence = T)

  outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)

  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop)
    return(outann)
  }
}


#' summary.acproj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the 'acproj()' function.
#' @param printpred Whether to print the observed and predicted number of cases (`T` or `F`).
#' @param printcall Whether to print function `Call` for acproj.object (`T` or `F`).
#' @param digits Number of digits in output, default (`0`) is integer only.
#' @param ... Other parameters.
#'
#' @return An information table describing the method used.
#'
#' @export
summary.acproj <- function(
  object,
  printpred = F,
  printcall = F,
  digits = 0,
  ...
) {
  method <- "Age-Cohort Model"

  if (!inherits(object, "acproj")) {
    stop("Variable \"object\" must be of type \"acproj\"")
  }

  obsto <- names(object$predictions)[
    dim(object$predictions)[2] - object$nopred
  ]

  if (!is.null(object$gofpvalue)) {
    gofpvalue <- round(object$gofpvalue, 4)
  } else {
    gofpvalue <- NA
  }

  if (printpred) {
    cat("Observed and predicted values:")
    cat("(observations up to", obsto, ")\n")
    print(round(as.matrix(object$predictions), digits = digits))
    cat("\n")
  }
  cat("\nPrediction done with:\n")

  moptions <- matrix(NA, 7, 2)
  moptions[, 1] <- c(
    "Method:",
    "Number of periods predicted (nopred):",
    "Trend used in new cohort estimation (cuttrend):",
    "Number of periods used in estimate (noperiod):",
    "Distribution function of regression:",
    "P-value for goodness of fit:",
    "First age group estimated (startestage):"
  )
  moptions[, 2] <- c(
    method,
    object$nopred,
    paste(object$cuttrend, collapse = " , "),
    object$noperiod,
    object$distribution,
    gofpvalue,
    object$startestage
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


#' glm.acproj
#'
#' Summarize estimations from the final model.
#'
#' @inheritParams acproj.getpred
#'
#' @return A summary table from the `glm.object`.
#'
#' @export
glm.acproj <- function(acproj.object) {
  summary(acproj.object$glm)
}


#' plot.acproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' acproj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
#' @param x An object based on the 'acproj()' function.
#'
#' @export
plot.acproj <- function(
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
  if (!inherits(x, "acproj")) {
    stop("Variable \"x\" must be of type \"acproj\"")
  }

  S7::check_is_S7(standpop, StandardPopulation)

  # Reading & formatting data:
  indat <- acproj.getproj(
    cdat,
    pdat,
    startp = startp,
    x,
    standpop = standpop
  )
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Setting internal variables:
  obsy <- dim(cdat)[2]
  nopredy <- length(indata) - obsy
  if (is.null(labels)) {
    labels <- row.names(indat)
  }

  # Create plots:
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

  # Returning object as invisible
  invisible(x)
}
