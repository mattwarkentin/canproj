#' ADPCPROJ
#'
#' R functions for projection of cancer incidence/mortality. Revising nordpred
#'   and introducing negative binomial distribution when lack of fit appears
#'   from nordpred, additional link functions of sqrt and identity, and settings
#'    of startestage and startuseage.
#'
#' @inheritParams canproj
#' @param n5case Minimum number of cancer cases/deaths per 5 years for splitting data.
#' @param noperiods List of candidate periods for projection base. Default (`NULL`)
#'  uses a goodness-of-fit test to determine if ancient periods are removed.
#' @param recent Estimate drift term from recent trend (`T`) or whole trend (`F`).
#'  Default (`NULL`) uses compares models to pick.
#'
#' @return A `list()`.
#'
#' @export
adpcproj <- function(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  noperiods = NULL,
  recent = NULL,
  startestage = NULL,
  newcohort = NULL,
  pGOF = 0.05,
  cuttrd = 0.04,
  shortp = 0,
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

  if (is.null(newcohort)) {
    startuseage <- startestage
  } else {
    startuseage <- startestage + 1
  }

  percases <- ncol(cases)
  if (percases < 3) {
    stop("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  if (is.null(noperiods)) {
    noperiods <- c(min(percases, 4):min(percases, 6))
  }

  noperiods <- sort(noperiods)
  while (length(noperiods) > 1) {
    maxnoperiod <- max(noperiods)
    glmn <- adpcproj.estimate(cases, pyr, maxnoperiod, startestage)$glm
    pval <- 1 - stats::pchisq(glmn$deviance, glmn$df.residual)
    if (pval < 0.01) {
      noperiods <- noperiods[1:(length(noperiods) - 1)]
    } else {
      noperiods <- maxnoperiod
    }
  }
  noperiod <- noperiods

  if (is.null(recent)) {
    recent <- adpcproj.estimate(
      cases,
      pyr,
      noperiod,
      startestage
    )$suggestionrecent
  }

  est <- adpcproj.estimate(
    cases = cases,
    pyr = pyr,
    noperiod = noperiod,
    pGOF = pGOF,
    startestage = startestage,
    linkfunc = linkfunc
  )
  pred <- adpcproj.prediction(
    adpcproj.estimate.object = est,
    startuseage = startuseage,
    recent = recent,
    cuttrd = cuttrd,
    shortp = shortp
  )
  return(pred)
}


#' adpcproj.estimate
#'
#' Fit age-drift-period-cohort models
#'
#' @inheritParams canproj
#' @param cases `data.frame` with number of cases in `nagg`-year period by ascending age groups in row.
#' @param pyr `data.frame` with observed and projected population size in `nagg`-year.
#' @param noperiod Number of 5-year periods in historical data.
#'
#' @return A `list()`.
#'
#' @export
adpcproj.estimate <- function(
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

  apcdata <- apcdata[apcdata$Age >= startestage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]

  options(contrasts = c("contr.treatment", "contr.poly"))

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
      Cases ~
        as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
      data = apcdata,
      family = power5link
    )
  } else if (linkfunc == "log") {
    res.glm <- stats::glm(
      Cases ~
        as.factor(Age) +
          Period +
          as.factor(Period) +
          as.factor(Cohort) +
          offset(log(y)) -
          1,
      data = apcdata,
      family = stats::poisson(link = log)
    )
  } else if (linkfunc == "sqrt") {
    suppressWarnings(
      res.glm <- stats::glm(
        Cases / y ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
        data = apcdata,
        family = stats::poisson(link = sqrt)
      )
    )
  } else if (linkfunc == "identity") {
    suppressWarnings(
      res.glm <- stats::glm(
        Cases / y ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
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

      glmnb <- suppressWarnings(MASS::glm.nb(
        Cases ~
          as.factor(Age) +
            Period +
            as.factor(Period) +
            as.factor(Cohort) -
            1 +
            offset(log(y)),
        data = apcdata,
        link = log
      ))
      theta <- as.numeric(suppressWarnings(MASS::theta.md(
        apcdata$Cases,
        stats::fitted(glmnb),
        dfr = stats::df.residual(glmnb)
      )))
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
        Cases ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
        data = apcdata,
        family = nbpower5link
      )
    } else if (linkfunc == "log") {
      res.glm <- suppressWarnings(MASS::glm.nb(
        Cases ~
          as.factor(Age) +
            Period +
            as.factor(Period) +
            as.factor(Cohort) -
            1 +
            offset(log(y)),
        data = apcdata,
        link = log
      ))
    } else if (linkfunc == "sqrt") {
      res.glm <- MASS::glm.nb(
        Cases / y ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
        data = apcdata,
        link = sqrt
      )
    } else if (linkfunc == "identity") {
      res.glm <- MASS::glm.nb(
        Cases / y ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
        data = apcdata,
        link = identity
      )
    } else {
      stop("Unknown \"linkfunc\"")
    }

    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  }

  if (distn == "Poisson") {
    mod1 <- stats::glm(
      Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 1,
      data = apcdata,
      family = stats::poisson
    )
    mod2 <- stats::glm(
      Cases ~
        as.factor(Age) +
          Period +
          I(Period^2) +
          as.factor(Cohort) +
          offset(log(y)) -
          1,
      data = apcdata,
      family = stats::poisson
    )
    pdiff <- 1 -
      stats::pchisq(
        (mod1$deviance - mod2$deviance),
        (mod1$df.residual - mod2$df.residual)
      )

    if (pdiff < 0.05) {
      suggestionrecent <- T
    } else {
      suggestionrecent <- F
    }
  } else {
    mod1 <- suppressWarnings(MASS::glm.nb(
      Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 1,
      data = apcdata,
      link = log
    ))
    mod2 <- suppressWarnings(MASS::glm.nb(
      Cases ~
        as.factor(Age) +
          Period +
          I(Period^2) +
          as.factor(Cohort) +
          offset(log(y)) -
          1,
      data = apcdata,
      link = log
    ))
    pdiff <- 1 -
      stats::pchisq(
        (mod2$twologlik - mod1$twologlik),
        (mod1$df.residual - mod2$df.residual)
      )
    if (pdiff < 0.05) {
      suggestionrecent <- T
    } else {
      suggestionrecent <- F
    }
  }

  res <- list(
    glm = res.glm,
    cases = cases,
    pyr = pyr,
    noperiod = noperiod,
    linkfunc = linkfunc,
    startestage = startestage,
    distribution = distn,
    gofpvalue = pvalue,
    suggestionrecent = suggestionrecent,
    pvaluerecent = pdiff
  )
  class(res) <- "adpcproj.estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' adpcproj.prediction
#'
#' Extrapolate estimated trend from age-drift-period-cohort model
#'
#' @inheritParams canproj
#' @param adpcproj.estimate.object An object based on the `adpcproj.estimate()` function.
#' @param startuseage Youngest age group to use estimates from the GLM for projection.
#' @param recent Indicate estimated drift term from recent trend (`T`) or whole trend (`F`).
#'
#' @return A `list()`.
#'
#' @export
adpcproj.prediction <- function(
  adpcproj.estimate.object,
  startuseage,
  recent,
  shortp = 0,
  cuttrd = 0.04
) {
  if (!inherits(adpcproj.estimate.object, "adpcproj.estimate")) {
    stop(
      "Variable \"adpcproj.estimate.object\" must be of type \"adpcproj.estimate\""
    )
  }

  if (adpcproj.estimate.object$startestage > startuseage) {
    stop("\"startuseage\" is set too low compared to \"startestage\"")
  }

  cases <- adpcproj.estimate.object$cases
  pyr <- adpcproj.estimate.object$pyr
  noperiod <- adpcproj.estimate.object$noperiod
  nototper <- ncol(pyr)
  noobsper <- ncol(cases)
  nonewpred <- nototper - noobsper
  cuttrend <- rep(shortp, nonewpred)
  ngroups <- nrow(cases)

  for (i in 1:nonewpred) {
    cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
  }
  cuttrend[cuttrend > 1] <- 1

  datatable <- data.frame(cases, matrix(NA, ngroups, nonewpred))
  names(datatable) <- names(pyr)

  skipped_ages <- 1:(startuseage - 1)
  obsinc <- cases[skipped_ages, (noobsper - 1):noobsper] /
    pyr[skipped_ages, (noobsper - 1):noobsper]
  obsinc[is.na(obsinc)] <- 0

  datatable[skipped_ages, (noobsper + 1):nototper] <- ((obsinc[, 1] +
    obsinc[, 2]) /
    2) *
    pyr[skipped_ages, (noobsper + 1):nototper]

  # For old agegroups, use trend from model estimates:
  for (age in startuseage:ngroups) {
    startestage <- adpcproj.estimate.object$startestage
    coefficients <- adpcproj.estimate.object$glm$coefficients

    coh <- (ngroups - startestage) -
      (age - startestage) +
      (noperiod + 1:nonewpred)

    noages <- ngroups - startestage + 1
    driftmp <- cumsum(1 - cuttrend)
    cohfind <- noages + (noperiod - 1) + 1 + (coh - 1)
    maxcoh <- ngroups - startuseage + noperiod
    agepar <- as.numeric(coefficients[age - startestage + 1])
    driftfind <- pmatch("Period", attributes(coefficients)$names)
    driftpar <- as.numeric(coefficients[driftfind])

    cohpar <- rep(NA, length(coh))
    for (i in 1:length(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        cohpar[i] <- as.numeric(coefficients[
          length(coefficients) - (startuseage - startestage)
        ])
        cohpar[i][is.na(cohpar[i])] <- 0
      }
    }

    if (recent) {
      lpfind <- driftfind + noperiod - 2
      lppar <- as.numeric(coefficients[lpfind])
      driftrecent <- driftpar - lppar
    }

    if (adpcproj.estimate.object$linkfunc == "power5") {
      if (recent) {
        rate <- (agepar +
          driftpar * noobsper +
          driftrecent * driftmp +
          cohpar)^5
      } else {
        rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)^5
      }
    } else if (adpcproj.estimate.object$linkfunc == "log") {
      if (recent) {
        rate <- exp(
          agepar + driftpar * noobsper + driftrecent * driftmp + cohpar
        )
      } else {
        rate <- exp(agepar + driftpar * (noobsper + driftmp) + cohpar)
      }
    } else if (adpcproj.estimate.object$linkfunc == "sqrt") {
      if (recent) {
        rate <- (agepar +
          driftpar * noobsper +
          driftrecent * driftmp +
          cohpar)^2
      } else {
        rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)^2
      }
    } else {
      # identity link:
      if (recent) {
        rate <- (agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)
      } else {
        rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)
      }
    }
    datatable[age, (noobsper + 1):nototper] <- rate *
      pyr[age, (noobsper + 1):nototper]
  }

  res <- list(
    predictions = datatable,
    pyr = pyr,
    nopred = nonewpred,
    noperiod = adpcproj.estimate.object$noperiod,
    gofpvalue = adpcproj.estimate.object$gofpvalue,
    recent = recent,
    pvaluerecent = adpcproj.estimate.object$pvaluerecent,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    distribution = adpcproj.estimate.object$distribution,
    startuseage = startuseage,
    startestage = adpcproj.estimate.object$startestage,
    glm = adpcproj.estimate.object$glm
  )
  class(res) <- "adpcproj"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' adpcproj.getpred
#'
#' Extract projection results by 5-year period
#'
#' @inheritParams canproj
#' @param adpcproj.object An object based on the `adpcproj()` function.
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
adpcproj.getpred <- function(
  adpcproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
) {
  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  if (!inherits(adpcproj.object, "adpcproj")) {
    stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")
  }

  validate_getpred_inputs(byage, standpop, incidence, agegroups)

  res <- make_pred_table(
    adpcproj.object,
    agegroups,
    standpop,
    byage,
    incidence,
    excludeobs
  )

  return(res)
}


#' adpcproj.getproj
#'
#' Extract annual projection results
#'
#' @inheritParams canproj
#' @inheritParams adpcproj.getpred
#'
#' @return A `data.frame()`.
#'
#' @export
adpcproj.getproj <- function(
  cdat,
  pdat,
  startp,
  adpcproj.object,
  standpop = NULL
) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  r0 <- adpcproj.getpred(adpcproj.object, incidence = T)
  outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)
  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop)
    return(outann)
  }
}


#' summary.adpcproj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the 'adpcproj()' function.
#' @param printpred Whether to print the observed and predicted number of cases (`T` or `F`).
#' @param printcall Whether to print function `Call` for adpcproj.object (`T` or `F`).
#' @param digits Number of digits in output, default (`0`) is integer only.
#' @param ... Other parameters.
#'
#' @return An information table describing the method used.
#'
#' @export
summary.adpcproj <- function(
  object,
  printpred = F,
  printcall = F,
  digits = 0,
  ...
) {
  method <- "Age-drift-Period-Cohort Model"
  if (!inherits(object, "adpcproj")) {
    stop("Variable \"object\" must be of type \"adpcproj\"")
  }
  obsto <- names(object$predictions)[
    ncol(object$predictions) - object$nopred
  ]

  if (!is.null(object$pvaluerecent)) {
    precent <- round(object$pvaluerecent, 4)
  } else {
    precent <- NA
  }

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

  moptions <- matrix(NA, 10, 2)
  moptions[, 1] <- c(
    "Method:",
    "Number of periods predicted (nopred):",
    "Trend used in predictions (cuttrend):",
    "Number of periods used in estimate (noperiod):",
    "Distribution function of regression:",
    "P-value for goodness of fit:",
    "Used recent (recent):",
    "P-value for recent:",
    "First age group used (startuseage):",
    "First age group estimated (startestage):"
  )
  moptions[, 2] <- c(
    method,
    object$nopred,
    paste(object$cuttrend, collapse = " , "),
    object$noperiod,
    object$distribution,
    gofpvalue,
    object$recent,
    precent,
    object$startuseage,
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


#' glm.adpcproj
#'
#' Summarize estimations from the final model.
#'
#' @inheritParams adpcproj.getpred
#'
#' @return A summary table from the `glm.object`.
#'
#' @export
glm.adpcproj <- function(adpcproj.object) {
  summary(adpcproj.object$glm)
}


#' plot.adpcproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' adpcproj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
#' @param x An object based on the 'adpcproj()' function.
#'
#' @export
plot.adpcproj <- function(
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
  if (!inherits(x, "adpcproj")) {
    stop("Variable \"x\" must be of type \"adpcproj\"")
  }

  S7::check_is_S7(standpop, StandardPopulation)

  # Reading & formatting data:
  indat <- adpcproj.getproj(
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
