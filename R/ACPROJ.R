#' ACPROJ
#'
#' R functions for projection of cancer incidence/mortality. Revising and
#'   combining nordpred and Osmond's to extrapolation cohort when no drift
#'   appears from nordpred.
#'
#' @inheritParams canproj
#' @param n5case Minimum number of cancer cases/deaths per 5 years for splitting data.
#' @param noperiods List of candidate periods for projection base. Default (`NULL`)
#'  uses a goodness-of-fit test to determine if ancient periods are removed.
#'
#' @export
acproj <- function(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  noperiods = NULL,
  startestage = NULL,
  cuttrd = 0.04,
  shortp = 0,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  # Define number of cases for data splitting:
  if (is.null(n5case)) {
    if (projfor == "incidence") {
      n5case <- 5
    } else {
      n5case <- 3
    }
  }

  ## aggregating data by 5 years:
  aggdata <- datagg(cdat, pdat, 5)
  cases <- aggdata$cases
  pyr <- aggdata$pyr

  ## Setting startestage:
  if (is.null(startestage)) {
    dat <- as.matrix(cases)
    iage <- 1
    while (mean(dat[iage, ]) < n5case | dat[iage, 1] == 0) {
      iage <- iage + 1
    }
    startestage <- iage
  }

  ## Setting default and checking data:
  # Number of periods in observed data
  percases <- dim(cases)[2]
  if (percases < 3) {
    stop("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  ## Setting number of periods for projection base
  noperiod <- percases

  ## Perform estimation and prediction:
  est <- acproj.estimate(
    cases = cases,
    pyr = pyr,
    noperiod = noperiod,
    pGOF = pGOF,
    startestage = startestage,
    linkfunc = linkfunc
  )

  pred <- acproj.prediction(
    acproj.estimate.object = est,
    cuttrd = 0.04,
    shortp = 0
  )

  return(pred)
}


## Fitting models:
#' @rdname acproj
#' @export
acproj.estimate <- function(
  cases,
  pyr,
  noperiod,
  startestage,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  ## Checking data
  if (dim(cases)[1] != 19 || dim(pyr)[1] != 19) {
    stop("\"cases\" and \"pyr\" must have data for 19 age groups")
  }

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

  ## Setting internal variables:
  dnoperiods <- dim(cases)[2]
  dnoagegr <- dim(cases)[1]

  ## Transform dataformat:
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

  ## Estimation:
  if (linkfunc == "power5") {
    ## Creation of power5 link:
    # Setting population variable
    y <- apcdata$y
    # Make power5 link function for poisson family:
    power5link <- poisson()
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
    res.glm <- glm(
      Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
      data = apcdata,
      family = power5link
    )
  } else if (linkfunc == "log") {
    res.glm <- glm(
      Cases ~ factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
      data = apcdata,
      family = poisson(link = log)
    )
  } else if (linkfunc == "sqrt") {
    res.glm <- glm(
      Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
      data = apcdata,
      family = poisson(link = sqrt)
    )
  } else if (linkfunc == "identity") {
    res.glm <- glm(
      Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
      data = apcdata,
      family = poisson(link = identity)
    )
  } else {
    stop("Unknown \"linkfunc\"")
  }
  pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

  ## change model to negative binomial glm when lack of fit:
  if (pvalue < pGOF) {
    distn <- "Negative-binomial"
    options(warn = -1)
    # Estimation:
    if (linkfunc == "power5") {
      y <- apcdata$y
      # Make power5 link function for negative binomial family:
      glmnb <- MASS::glm.nb(
        Cases ~
          factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
        data = apcdata,
        link = log
      )
      theta <- as.numeric(MASS::theta.md(
        apcdata$Cases,
        fitted(glmnb),
        dfr = df.residual(glmnb)
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
      res.glm <- glm(
        Cases ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        family = nbpower5link
      )
    } else if (linkfunc == "log") {
      res.glm <- MASS::glm.nb(
        Cases ~
          factor(Age) + relevel(factor(Cohort), midc) + offset(log(y)) - 1,
        data = apcdata,
        link = log
      )
    } else if (linkfunc == "sqrt") {
      res.glm <- MASS::glm.nb(
        Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        link = sqrt
      )
    } else if (linkfunc == "identity") {
      res.glm <- MASS::glm.nb(
        Cases / y ~ factor(Age) + relevel(factor(Cohort), midc) - 1,
        data = apcdata,
        link = identity
      )
    } else {
      stop("Unknown \"linkfunc\"")
    }
    # updating p-value of goodness of fit and distribution function:
    pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
    options(warn = 0)
  }

  ## Set class and return results
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


## Project age-specific rates:
#' @rdname acproj
#' @export
acproj.prediction <- function(
  acproj.estimate.object,
  cuttrd = 0.04,
  shortp = 0
) {
  ## running conditions:
  if (class(acproj.estimate.object) != "acproj.estimate") {
    stop(
      "Variable \"acproj.estimate.object\" must be of type \"acproj.estimate\""
    )
  }

  ## Setting local variables:
  cases <- acproj.estimate.object$cases
  pyr <- acproj.estimate.object$pyr
  noperiod <- acproj.estimate.object$noperiod
  startestage <- acproj.estimate.object$startestage
  maxc <- acproj.estimate.object$maxc
  midc <- acproj.estimate.object$midc
  nototper <- dim(pyr)[2]
  noobsper <- dim(cases)[2]
  nonewpred <- nototper - noobsper

  cuttrend <- rep(shortp, nonewpred)
  for (i in 1:nonewpred) {
    cuttrend[i] <- shortp + (i - 1) * 5 * cuttrd
  }
  cuttrend[cuttrend > 1] <- 1
  if (length(cuttrend) > nonewpred) {
    cuttrend <- cuttrend[1:nonewpred]
  }

  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode", 1:nototper)
    } else {
      years <- dimnames(pyr)[[2]]
    }
  }

  ## Making data object:
  datatable <- matrix(NA, 19, nototper)
  # fill in observed cases:
  datatable[, 1:(nototper - nonewpred)] <- as.matrix(cases)
  datatable <- data.frame(datatable)
  row.names(datatable) <- c(
    "0-4",
    "5-9",
    "10-14",
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75-79",
    "80-84",
    "85-89",
    "90+"
  )
  names(datatable) <- years

  ## Calculate predictions in number of cases:
  # For young agegroups with little data, use average from last two periods:
  for (age in 1:(startestage - 1)) {
    obsinc <- cases[age, (noobsper - 1):noobsper] /
      pyr[age, (noobsper - 1):noobsper]
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0
    }
    datatable[age, (noobsper + 1):nototper] <- ((obsinc[, 1] + obsinc[, 2]) /
      2) *
      pyr[age, (noobsper + 1):nototper]
  }

  # For old agegroups, use AC model estimates:
  coefficients <- acproj.estimate.object$glm$coefficients
  # Age effects:
  noages <- 19 - startestage + 1
  ageff <- coefficients[1:noages]
  # Cohort effects:
  coeff <- c(
    coefficients[(noages + 1):(noages + midc - 1)],
    0,
    coefficients[(noages + midc):length(coefficients)]
  )
  lncoh <- length(coeff)
  # For new young cohorts in "nonewpred" by linear regression without the first and last 4 cohorts:
  ceff <- coeff[-c(1, 2, lncoh - 1, lncoh)]
  lnceff <- length(ceff)
  tt <- 1:lnceff
  coeflm <- coef(lm(ceff ~ tt - 1))
  ncoeff <- coeflm *
    c((lncoh - 3):(lncoh - 2), (cumsum(1 - cuttrend) + lncoh - 2))
  coheff <- c(coeff[1:(lncoh - 2)], ncoeff) # replace the last 2 cohort by estimates

  for (age in startestage:19) {
    # Age effect:
    agepar <- as.numeric(ageff[age - startestage + 1])
    # Cohort index: No. agegoups - age + period
    coh <- (19 - startestage) - (age - startestage) + (noperiod + 1:nonewpred)
    cohpar <- coheff[coh]
    # Project age-specific rates:
    if (acproj.estimate.object$linkfunc == "power5") {
      rate <- (agepar + cohpar)^5
    } else if (acproj.estimate.object$linkfunc == "log") {
      rate <- exp(agepar + cohpar)
    } else if (acproj.estimate.object$linkfunc == "sqrt") {
      rate <- (agepar + cohpar)^2
    } else {
      # identity link:
      rate <- (agepar + cohpar)
    }
    datatable[age, (noobsper + 1):nototper] <- rate *
      pyr[age, (noobsper + 1):nototper]
  }

  # Structure and return results:
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


## Summary Projection Results by 5-year Period(ASR and total number):
#' @rdname acproj
#' @export
acproj.getpred <- function(
  acproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
) {
  ## Setting defaults:
  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  ## Checking input:
  if (class(acproj.object) != "acproj") {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")
  }

  if ((!is.null(standpop)) && (!incidence)) {
    stop(
      "\"standpop\" should only be used with incidence predictions (incidence=T)"
    )
  }

  if (!is.null(standpop)) {
    if (round(sum(standpop), 5) != 1) {
      stop("\"standpop\" must be of sum 1")
    }
    if ((length(standpop) != length(agegroups)) && (agegroups[1] != "all")) {
      stop("\"standpop\" must be the same length as \"agegroups\"")
    }
    if (byage) {
      stop("\"standpop\" is only valid for \"byage=T\"")
    }
  }

  ## Seting local data:
  datatable <- acproj.object$predictions
  pyr <- data.frame(acproj.object$pyr)

  ## Secting agegroups:
  if (agegroups[1] != "all") {
    datatable <- datatable[agegroups, ]
    pyr <- pyr[agegroups, ]
  }

  ## If needed; Standardize data and Collapse agegroups
  if (!is.null(standpop)) {
    datainc <- (datatable / pyr) * 100000
    if (sum(is.na(datainc)) > 0) {
      datainc[is.na(datainc)] <- 0
    }
    res <- apply(datainc * standpop, 2, sum)
  } else {
    if (!byage) {
      datatable <- apply(datatable, 2, sum)
      pyr <- apply(pyr, 2, sum)
    }
    if (incidence) {
      res <- (datatable / pyr) * 100000
      if (sum(is.na(res)) > 0) {
        res[is.na(res)] <- 0
      }
    } else {
      res <- datatable
    }
  }

  ## Select data:
  if (excludeobs) {
    if (is.matrix(res)) {
      predstart <- dim(res)[2] - acproj.object$nopred + 1
      res <- res[, predstart:(predstart + acproj.object$nopred - 1)]
    } else {
      predstart <- length(res) - acproj.object$nopred + 1
      res <- res[predstart:(predstart + acproj.object$nopred - 1)]
    }
  }

  ## Return data:
  return(res)
}


## Summary Annual Projection Results:
#' @rdname acproj
#' @export
acproj.getproj <- function(cdat, pdat, startp, acproj.object, standpop = NULL) {
  r0 <- acproj.getpred(acproj.object, incidence = T)

  outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)

  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop = standpop)
    return(outann)
  }
}


## Summaring methods used for the projection:
#' @rdname acproj
#' @export
summary.acproj <- function(
  acproj.object,
  printpred = F,
  printcall = F,
  digits = 0
) {
  method <- "Age-Cohort Model"

  if (class(acproj.object) != "acproj") {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")
  }

  # Setting internal variables:
  obsto <- names(acproj.object$predictions)[
    dim(acproj.object$predictions)[2] - acproj.object$nopred
  ]

  if (!is.null(acproj.object$gofpvalue)) {
    gofpvalue <- round(acproj.object$gofpvalue, 4)
  } else {
    gofpvalue <- NA
  }

  # Print information about object:
  if (printpred) {
    cat("Observed and predicted values:")
    cat("(observations up to", obsto, ")\n")
    print(round(as.matrix(acproj.object$predictions), digits = digits))
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
    acproj.object$nopred,
    paste(acproj.object$cuttrend, collapse = " , "),
    acproj.object$noperiod,
    acproj.object$distribution,
    gofpvalue,
    acproj.object$startestage
  )

  maxl <- max(nchar(moptions[, 1]))

  for (i in 1:dim(moptions)[1]) {
    spaces <- rep(" ", maxl - nchar(moptions[i, 1]) + 2)
    cat(moptions[i, 1], spaces, moptions[i, 2], "\n", sep = "")
  }
  if (printcall) {
    cat("\n  Call: ")
    dput(attr(acproj.object, "Call"))
  }
  invisible(acproj.object)
}

## Obtaining modeling info and parameter estimates:
#' @rdname acproj
#' @export
glm.acproj <- function(acproj.object) {
  summary(acproj.object$glm)
}

## Ploting prejected standardized rates:
#' @rdname acproj
#' @export
plot.acproj <- function(
  cdat,
  pdat,
  startp,
  acproj.object,
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
  if (class(acproj.object) != "acproj") {
    stop("Variable \"acproj.object\" must be of type \"acproj\"")
  }

  # Reading & formating data:
  indat <- acproj.getproj(
    cdat,
    pdat,
    startp = startp,
    acproj.object,
    standpop = standpop
  )
  indata <- indat[, 1]
  indata <- indata[startplot:length(indata)]

  # Seting internal variables:
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
    plot(c(1, maxx), ylim, type = "n", ylab = ylab, xlab = xlab, axes = F, ...)
    axis(2)
    axis(1, at = 1:maxx, labels = labels)
    box()
    title(main)
  }

  lines(
    1:(maxx - nopredy),
    indata[1:(maxx - nopredy)],
    type = "o",
    pch = 20,
    lty = lty[1],
    col = col[1],
    ...
  )

  lines(
    (maxx - nopredy):maxx,
    indata[(maxx - nopredy):maxx],
    lty = lty[2],
    col = col[2],
    ...
  )

  # Returning object as invisible
  invisible(acproj.object)
}
