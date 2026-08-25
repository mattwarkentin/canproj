#' ACPROJ: Age-Cohort Projections
#'
#' R functions for projection of cancer incidence/mortality. Revising and
#'   combining nordpred and Osmond's to extrapolation cohort when no drift
#'   appears from nordpred.
#'
#' @inheritParams canproj
#' @param n5case Minimum number of cancer cases/deaths per 5 years for splitting data.
#'
#' @return `acproj()` returns a `list`.
#'
#' @export
acproj <- function(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  startage = NULL,
  cuttrd = 0.04,
  shortp = 0,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  if (ncol(cdat) < 15) {
    rlang::abort("\"cdat\" must have at least 15 years")
  }

  if (ncol(pdat) < 20) {
    rlang::abort("\"pdat\" must have at least 20 years")
  }

  validate_projection_inputs(
    cdat,
    pdat,
    projfor = projfor,
    n5case = n5case,
    linkfunc = linkfunc,
    startage = startage,
    cuttrd = cuttrd,
    shortp = shortp,
    pGOF = pGOF
  )

  if (is.null(n5case)) {
    if (projfor == "incidence") {
      n5case <- 5
    } else {
      n5case <- 3
    }
  }

  aggdata <- aggregate_data(cdat, pdat, 5)
  cases <- aggdata$cases
  pyr <- aggdata$pyr

  if (is.null(startage)) {
    dat <- as.matrix(cases)
    iage <- 1
    while (mean(dat[iage, ]) < n5case | dat[iage, 1] == 0) {
      iage <- iage + 1
    }
    startage <- iage
  }

  percases <- ncol(cases)
  if (percases < 3) {
    rlang::abort("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  est <- acproj_estimate(
    cases = cases,
    pyr = pyr,
    noperiod = percases,
    pGOF = pGOF,
    startage = startage,
    linkfunc = linkfunc
  )

  pred <- acproj_predict(
    acproj_estimate.object = est,
    cuttrd = cuttrd,
    shortp = shortp
  )

  return(pred)
}


#' acproj_estimate
#'
#' Fit age-cohort models
#'
#' @inheritParams canproj
#' @param cases `data.frame` with number of cases in `nagg`-year period by ascending age groups in row.
#' @param pyr `data.frame` with observed and projected population size in `nagg`-year.
#' @param noperiod Number of 5-year periods in historical data.
#'
#' @return A `list()`.
#' @keywords internal
acproj_estimate <- function(
  cases,
  pyr,
  noperiod,
  startage,
  pGOF = 0.05,
  linkfunc = "power5"
) {
  if ((ncol(pyr) - ncol(cases)) > 6) {
    rlang::abort("Package can not project more than 6 periods")
  }

  if (!inherits(noperiod, "numeric") & !inherits(noperiod, "integer")) {
    rlang::abort(
      "Variable \"noperiod\" must be of type \"numeric\" or \"integer\""
    )
  }

  if ((ncol(cases) - noperiod) < 0) {
    rlang::abort(
      "More periods specified in \"noperiod\" than columns in \"cases\""
    )
  }

  if (noperiod < 3) {
    rlang::abort("\"noperiod\" must be 3 or larger")
  }

  if (!inherits(startage, "numeric")) {
    rlang::abort("Variable \"startage\" must be of type \"numeric\"")
  }

  if (startage < 1) {
    rlang::abort("Variable \"startage\" must be >= 1")
  }

  if (startage > nrow(cases) - 1) {
    rlang::abort(
      "Variable \"startage\" is too high relative to number of age groups"
    )
  }

  validate_estimate_inputs(pyr, cases, pGOF, linkfunc)

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
  apcdata <- apcdata[apcdata$Age >= startage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]
  maxc <- max(apcdata$Cohort)
  midc <- ceiling(maxc / 2)

  formula_string <- paste0(
    "Cases ~ factor(Age) + relevel(factor(Cohort), ",
    midc,
    ") - 1"
  )
  res.glm <- get_glm(formula_string, apcdata, linkfunc)

  pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

  if (pvalue < pGOF) {
    distn <- "Negative-binomial"

    formula_string <- paste0(
      "Cases ~ factor(Age) + relevel(factor(Cohort), ",
      midc,
      ") - 1"
    )
    res.glm <- get_glm_nb(formula_string, apcdata, linkfunc)

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
    startage = startage,
    distribution = distn,
    gofpvalue = pvalue
  )
  class(res) <- "acproj_estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' acproj_predict
#'
#' Extrapolate estimated trend from age-cohort model
#'
#' @inheritParams canproj
#' @param acproj_estimate.object An object based on the `acproj_estimate()` function.
#'
#' @return A `list()`.
#' @keywords internal
acproj_predict <- function(
  acproj_estimate.object,
  cuttrd = 0.04,
  shortp = 0
) {
  if (!inherits(acproj_estimate.object, "acproj_estimate")) {
    rlang::abort(
      "Variable \"acproj_estimate.object\" must be of type \"acproj_estimate\""
    )
  }

  validate_prediction_inputs(cuttrd, shortp)

  cases <- acproj_estimate.object$cases
  pyr <- acproj_estimate.object$pyr
  noperiod <- acproj_estimate.object$noperiod
  startage <- acproj_estimate.object$startage
  maxc <- acproj_estimate.object$maxc
  midc <- acproj_estimate.object$midc
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

  skipped_ages <- 1:(startage - 1)
  obsinc <- cases[skipped_ages, (noobsper - 1):noobsper] /
    pyr[skipped_ages, (noobsper - 1):noobsper]
  obsinc[is.na(obsinc)] <- 0

  datatable[skipped_ages, (noobsper + 1):nototper] <- ((obsinc[, 1] +
    obsinc[, 2]) /
    2) *
    pyr[skipped_ages, (noobsper + 1):nototper]

  coefficients <- acproj_estimate.object$glm$coefficients

  num_ages <- ngroups - startage + 1
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

  for (age in startage:ngroups) {
    agepar <- as.numeric(ageff[age - startage + 1])
    coh <- (ngroups - startage) -
      (age - startage) +
      (noperiod + 1:nonewpred)
    cohpar <- coheff[coh]

    if (acproj_estimate.object$linkfunc == "power5") {
      rate <- (agepar + cohpar)^5
    } else if (acproj_estimate.object$linkfunc == "log") {
      rate <- exp(agepar + cohpar)
    } else if (acproj_estimate.object$linkfunc == "sqrt") {
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
    noperiod = acproj_estimate.object$noperiod,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    gofpvalue = acproj_estimate.object$gofpvalue,
    distribution = acproj_estimate.object$distribution,
    startage = acproj_estimate.object$startage,
    glm = acproj_estimate.object$glm
  )

  class(res) <- c("acproj", "proj_model")
  attr(res, "Call") <- sys.call()
  res
}


#' acproj_get_predictions
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
#' @keywords internal
acproj_get_predictions <- function(
  acproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage = NULL,
  agegroups = "all"
) {
  if (is.null(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  if (!inherits(acproj.object, "acproj")) {
    rlang::abort("Variable \"acproj.object\" must be of type \"acproj\"")
  }

  validate_getpred_inputs(byage, standpop, incidence, agegroups, excludeobs)

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


#' Get ACPROJ Projections
#'
#' `get_projections()` extracts annual projection results from an `acproj()`
#'   object.
#'
#' @inheritParams canproj
#' @inheritParams acproj_get_predictions
#' @param object An output object from `acproj()`.
#'
#' @return `get_projections()` returns a `data.frame`.
#'
#' @rdname acproj
#'
#' @export
get_projections.acproj <- function(
  object,
  ...,
  cdat,
  pdat,
  startp,
  standpop = NULL
) {
  rlang::check_dots_empty()
  validate_getproj_inputs(
    object = object,
    cdat = cdat,
    pdat = pdat,
    startp = startp,
    standpop = standpop
  )

  if (floor(ncol(cdat) / 5) != object$noperiod) {
    rlang::abort(
      "\"cdat\" must match acproj.object periods (floor(ncol(cdat) / 5) == noperiod)"
    )
  }

  if (floor(ncol(pdat) / 5) != ncol(object$predictions)) {
    rlang::abort(
      "\"pdat\" must match acproj.object periods (floor(ncol(pdat) / 5) == ncol(predictions))"
    )
  }

  r0 <- acproj_get_predictions(object, incidence = T)

  outasp <- interpolate_age_specific_rates(
    r0,
    cdat,
    pdat,
    startp = startp,
    nagg = 5
  )

  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- standardize_annual_rates(outasp, pdat, startp, standpop)
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
    rlang::abort("Variable \"object\" must be of type \"acproj\"")
  }

  validate_summary_inputs(
    printpred = printpred,
    printcall = printcall,
    digits = digits
  )

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
    "First age group estimated (startage):"
  )
  moptions[, 2] <- c(
    method,
    object$nopred,
    paste(object$cuttrend, collapse = " , "),
    object$noperiod,
    object$distribution,
    gofpvalue,
    object$startage
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


#' plot.acproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' acproj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
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
  ylab = "Rate per 100,000 people",
  main = "",
  lty = c(1, 3),
  col = c("black", "azure4"),
  ...
) {
  S7::check_is_S7(standpop, StandardPopulation)

  indat <- get_projections(
    object = x,
    cdat = cdat,
    pdat = pdat,
    startp = startp,
    standpop = standpop
  )
  indata <- indat[, 1]
  obsy <- 5 * x$noperiod

  data <- as.data.frame(indata)
  data$year <- as.numeric(names(indata))

  data$Period <- "Observed"
  data$Period[(obsy + 1):length(indata)] <- "Projected"
  dupe <- data[obsy, 1:2]
  dupe$Period <- "Projected"

  data <- data[startplot:nrow(data), ]
  data <- rbind(dupe, data)

  custom_colours <- c("Observed" = col[1], "Projected" = col[2])
  custom_line <- c("Observed" = lty[1], "Projected" = lty[2])

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$year, y = indata, color = .data$Period)
  ) +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$Period)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(values = custom_colours) +
    ggplot2::scale_linetype_manual(values = custom_line) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      plot.title = ggplot2::element_text(size = 15)
    )
}
