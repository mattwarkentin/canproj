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
  ## Setting startestage and startuseage:
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
    startuseage <- startestage + 1 # assign newcohort effects as the last estimated cohort effect
  }

  ## Setting default and checking data:
  # Number of periods in observed data
  percases <- dim(cases)[2]
  if (percases < 3) {
    stop("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  ## Setting number of periods for projection base
  # List possible candidates for number of periods to base predictions on
  # Default is 4:6 if available
  if (is.null(noperiods)) {
    noperiods <- c(min(percases, 4):min(percases, 6))
  }

  # Choose number of periods by cutting stepwise execution of the
  # highest candidate number (i.e. cutting the most ancient periods)
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

  ## Setting status for recent (whether to use recent trend or whole trend)
  if (is.null(recent)) {
    recent <- adpcproj.estimate(
      cases,
      pyr,
      noperiod,
      startestage
    )$suggestionrecent
  }

  ## Perform estimation and prediction:
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
    cuttrd = 0.04,
    shortp = 0
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

  ## Setting contrast for age+drift+period+cohort
  options(contrasts = c("contr.treatment", "contr.poly"))

  ## Estimation:
  if (linkfunc == "power5") {
    ## Creation of power5 link:
    # Setting population variable
    y <- apcdata$y
    # Make power5 link function for poisson family:
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
    res.glm <- stats::glm(
      Cases / y ~
        as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
      data = apcdata,
      family = stats::poisson(link = sqrt)
    )
  } else if (linkfunc == "identity") {
    res.glm <- stats::glm(
      Cases / y ~
        as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
      data = apcdata,
      family = stats::poisson(link = identity)
    )
  } else {
    stop("Unknown \"linkfunc\"")
  }
  pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
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
          as.factor(Age) +
            Period +
            as.factor(Period) +
            as.factor(Cohort) -
            1 +
            offset(log(y)),
        data = apcdata,
        link = log
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
        Cases ~
          as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1,
        data = apcdata,
        family = nbpower5link
      )
    } else if (linkfunc == "log") {
      res.glm <- MASS::glm.nb(
        Cases ~
          as.factor(Age) +
            Period +
            as.factor(Period) +
            as.factor(Cohort) -
            1 +
            offset(log(y)),
        data = apcdata,
        link = log
      )
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
    # updating p-value of goodness of fit and distribution function:
    pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
    options(warn = 0)
  }

  ## setting suggestion for 'recent' (whether to use recent trend or whole trend)
  if (distn == "poisson") {
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
    options(warn = -1)
    mod1 <- MASS::glm.nb(
      Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 1,
      data = apcdata,
      link = log
    )
    mod2 <- MASS::glm.nb(
      Cases ~
        as.factor(Age) +
          Period +
          I(Period^2) +
          as.factor(Cohort) +
          offset(log(y)) -
          1,
      data = apcdata,
      link = log
    )
    options(warn = 0)
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

  #options(contrasts("unordered", "ordered"))

  ## Set class and return results
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
  ## running conditions:
  if (!inherits(adpcproj.estimate.object, "adpcproj.estimate")) {
    stop(
      "Variable \"adpcproj.estimate.object\" must be of type \"adpcproj.estimate\""
    )
  }

  if (adpcproj.estimate.object$startestage > startuseage) {
    stop("\"startuseage\" is set too high compared to \"startestage\"")
  }

  ## Setting local variables:
  cases <- adpcproj.estimate.object$cases
  pyr <- adpcproj.estimate.object$pyr
  noperiod <- adpcproj.estimate.object$noperiod
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
  if (
    length(cuttrend) <
      (dim(adpcproj.estimate.object$pyr)[2] -
        dim(adpcproj.estimate.object$cases)[2])
  ) {
    err <- paste("\"cuttrend\" must always be at least the same length as")
    err <- paste(err, "the number of periods with population forecasts")
    stop(err)
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
  for (age in 1:(startuseage - 1)) {
    obsinc <- cases[age, (noobsper - 1):noobsper] /
      pyr[age, (noobsper - 1):noobsper]
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0
    }
    datatable[age, (noobsper + 1):nototper] <- ((obsinc[, 1] + obsinc[, 2]) /
      2) *
      pyr[age, (noobsper + 1):nototper]
  }

  # For old agegroups, use trend from model estimates:
  for (age in startuseage:19) {
    startestage <- adpcproj.estimate.object$startestage
    coefficients <- adpcproj.estimate.object$glm$coefficients

    # Cohort index: No. agegoups - age + period
    coh <- (19 - startestage) - (age - startestage) + (noperiod + 1:nonewpred)

    noages <- 19 - startestage + 1
    driftmp <- cumsum(1 - cuttrend)
    cohfind <- noages + (noperiod - 1) + 1 + (coh - 1)
    maxcoh <- 19 - startuseage + noperiod
    agepar <- as.numeric(coefficients[age - startestage + 1])
    driftfind <- pmatch("Period", attributes(coefficients)$names)
    driftpar <- as.numeric(coefficients[driftfind])

    cohpar <- rep(NA, length(coh))
    for (i in 1:length(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        # For new young cohorts in the future, define the effect as the last cohort:
        cohpar[i] <- as.numeric(coefficients[
          length(coefficients) - (startuseage - startestage)
        ])
        #       cohpar[i][is.na(cohpar[i])] <- as.numeric(coefficients[length(coefficients)-1])
        cohpar[i][is.na(cohpar[i])] <- 0
      }
    }

    # Getting recent drift (D_last) estimate:
    if (recent) {
      # localize the last period effect (note: p.first=p.last)
      lpfind <- driftfind + noperiod - 2
      # find the estimate of the last two period effect (P-1)
      lppar <- as.numeric(coefficients[lpfind])
      driftrecent <- driftpar - lppar
    }

    # Project age-specific rates:
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

  # Structure and return results:
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
  ## Setting defaults:
  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  ## Checking input:
  if (!inherits(adpcproj.object, "adpcproj")) {
    stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")
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
  datatable <- adpcproj.object$predictions
  pyr <- data.frame(adpcproj.object$pyr)

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
      predstart <- dim(res)[2] - adpcproj.object$nopred + 1
      res <- res[, predstart:(predstart + adpcproj.object$nopred - 1)]
    } else {
      predstart <- length(res) - adpcproj.object$nopred + 1
      res <- res[predstart:(predstart + adpcproj.object$nopred - 1)]
    }
  }

  ## Return data:
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
  r0 <- adpcproj.getpred(adpcproj.object, incidence = T)
  outasp <- asrpy(r0, cdat, pdat, startp = startp, nagg = 5)
  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- asry(outasp, pdat, standpop = standpop)
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
  # Setting internal variables:
  obsto <- names(object$predictions)[
    dim(object$predictions)[2] - object$nopred
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

  # Print information about object:
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
#' @inheritParams adpcproj.getpred
#'
#' @export
plot.adpcproj <- function(
  cdat,
  pdat,
  startp,
  adpcproj.object,
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
  if (!inherits(adpcproj.object, "adpcproj")) {
    stop("Variable \"adpcproj.object\" must be of type \"adpcproj\"")
  }

  # Reading & formating data:
  indat <- adpcproj.getproj(
    cdat,
    pdat,
    startp = startp,
    adpcproj.object,
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
  invisible(adpcproj.object)
}
