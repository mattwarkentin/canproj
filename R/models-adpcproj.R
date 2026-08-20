#' ADPCPROJ: Age-Drift-Period-Cohort Projections
#'
#' R functions for projection of cancer incidence/mortality. Revising nordpred
#'   and introducing negative binomial distribution when lack of fit appears
#'   from nordpred, additional link functions of sqrt and identity, and settings
#'    of startage and startuseage.
#'
#' @inheritParams canproj
#' @param n5case Minimum number of cancer cases/deaths per 5 years for splitting data.
#' @param noperiods List of candidate periods for projection base. Default (`NULL`)
#'  uses a goodness-of-fit test to determine if ancient periods are removed.
#' @param recent Estimate drift term from recent trend (`T`) or whole trend (`F`).
#'  Default (`NULL`) uses compares models to pick.
#'
#' @return `adpcproj()` returns a `list`.
#'
#' @export
adpcproj <- function(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  noperiods = NULL,
  recent = NULL,
  startage = NULL,
  newcohort = FALSE,
  pGOF = 0.05,
  cuttrd = 0.04,
  shortp = 0,
  linkfunc = "power5"
) {
  if (ncol(cdat) < 15) {
    rlang::abort("\"cdat\" must have at least 15 years")
  }

  if (ncol(pdat) < 20) {
    rlang::abort("\"pdat\" must have at least 20 years")
  }

  if (
    !is.null(noperiods) &
      !inherits(noperiods, "integer") &
      !inherits(noperiods, "numeric")
  ) {
    rlang::abort(
      "\"noperiods\" must be of type \"integer\", \"numeric\", o
    r NULL"
    )
  }

  if (!is.null(recent) && !inherits(recent, "logical")) {
    rlang::abort("\"recent\" must be of type \"logical\" or NULL")
  }

  validate_projection_inputs(
    cdat,
    pdat,
    projfor = projfor,
    n5case = n5case,
    startage = startage,
    newcohort = newcohort,
    pGOF = pGOF,
    cuttrd = cuttrd,
    shortp = shortp,
    linkfunc = linkfunc
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

  if (!newcohort) {
    startuseage <- startage
  } else {
    startuseage <- startage + 1
  }

  percases <- ncol(cases)
  if (percases < 3) {
    rlang::abort("Minimum number of period is 3 (15 years) in \"cases\"")
  }

  if (is.null(noperiods)) {
    noperiods <- c(min(percases, 4):min(percases, 6))
  }

  noperiods <- sort(noperiods)
  while (length(noperiods) > 1) {
    maxnoperiod <- max(noperiods)
    glmn <- adpcproj_estimate(cases, pyr, maxnoperiod, startage)$glm
    pval <- 1 - stats::pchisq(glmn$deviance, glmn$df.residual)
    if (pval < 0.01) {
      noperiods <- noperiods[1:(length(noperiods) - 1)]
    } else {
      noperiods <- maxnoperiod
    }
  }
  noperiod <- noperiods

  if (is.null(recent)) {
    recent <- adpcproj_estimate(
      cases,
      pyr,
      noperiod,
      startage
    )$suggestionrecent
  }

  est <- adpcproj_estimate(
    cases = cases,
    pyr = pyr,
    noperiod = noperiod,
    pGOF = pGOF,
    startage = startage,
    linkfunc = linkfunc
  )
  pred <- adpcproj_predict(
    adpcproj_estimate.object = est,
    startuseage = startuseage,
    recent = recent,
    cuttrd = cuttrd,
    shortp = shortp
  )
  return(pred)
}


#' adpcproj_estimate
#'
#' Fit age-drift-period-cohort models
#'
#' @inheritParams canproj
#' @param cases `data.frame` with number of cases in `nagg`-year period by ascending age groups in row.
#' @param pyr `data.frame` with observed and projected population size in `nagg`-year.
#' @param noperiod Number of 5-year periods in historical data.
#'
#' @return A `list()`.
#' @keywords internal
adpcproj_estimate <- function(
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
      "Variable \"noperiod\" must be of type \"numeri
    c\" or \"integer\""
    )
  }

  if ((ncol(cases) - noperiod) < 0) {
    rlang::abort(
      "More periods specified in \"noperi
    od\" than columns in \"cases\""
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
      "Variable \"startage\" is
     too high relative to number of age groups"
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

  apcdata <- apcdata[apcdata$Age >= startage, ]
  apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]

  options(contrasts = c("contr.treatment", "contr.poly"))

  formula_string <- "Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1"
  res.glm <- get_glm(formula_string, apcdata, linkfunc)

  pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  distn <- "Poisson"

  if (pvalue < pGOF) {
    distn <- "Negative-binomial"

    formula_string <- "Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 1"
    res.glm <- get_glm_nb(formula_string, apcdata, linkfunc)

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
    startage = startage,
    distribution = distn,
    gofpvalue = pvalue,
    suggestionrecent = suggestionrecent,
    pvaluerecent = pdiff
  )
  class(res) <- "adpcproj_estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' adpcproj_predict
#'
#' Extrapolate estimated trend from age-drift-period-cohort model
#'
#' @inheritParams canproj
#' @param adpcproj_estimate.object An object based on the `adpcproj_estimate()` function.
#' @param startuseage Youngest age group to use estimates from the GLM for projection.
#' @param recent Indicate estimated drift term from recent trend (`T`) or whole trend (`F`).
#'
#' @return A `list()`.
#' @keywords internal
adpcproj_predict <- function(
  adpcproj_estimate.object,
  startuseage,
  recent,
  shortp = 0,
  cuttrd = 0.04
) {
  if (!inherits(adpcproj_estimate.object, "adpcproj_estimate")) {
    rlang::abort(
      "Variable \"adpcproj_estimate.object\" must be of type \"adpcproj_estimate\""
    )
  }

  if (adpcproj_estimate.object$startage > startuseage) {
    rlang::abort("\"startuseage\" is set too low compared to \"startage\"")
  }

  if (!inherits(recent, "logical")) {
    rlang::abort("\"recent\" must be of type \"logical\"")
  }

  if (!inherits(startuseage, "numeric")) {
    rlang::abort("\"startuseage\" must be of type \"numeric\"")
  }

  validate_prediction_inputs(cuttrd, shortp)

  cases <- adpcproj_estimate.object$cases
  pyr <- adpcproj_estimate.object$pyr
  noperiod <- adpcproj_estimate.object$noperiod
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

  for (age in startuseage:ngroups) {
    startage <- adpcproj_estimate.object$startage
    coefficients <- adpcproj_estimate.object$glm$coefficients

    coh <- (ngroups - startage) -
      (age - startage) +
      (noperiod + 1:nonewpred)

    noages <- ngroups - startage + 1
    driftmp <- cumsum(1 - cuttrend)
    cohfind <- noages + (noperiod - 1) + 1 + (coh - 1)
    maxcoh <- ngroups - startuseage + noperiod
    agepar <- as.numeric(coefficients[age - startage + 1])
    driftfind <- pmatch("Period", attributes(coefficients)$names)
    driftpar <- as.numeric(coefficients[driftfind])

    cohpar <- rep(NA, length(coh))
    for (i in 1:length(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        cohpar[i] <- as.numeric(coefficients[
          length(coefficients) - (startuseage - startage)
        ])
        cohpar[i][is.na(cohpar[i])] <- 0
      }
    }

    if (recent) {
      lpfind <- driftfind + noperiod - 2
      lppar <- as.numeric(coefficients[lpfind])
      driftrecent <- driftpar - lppar
    }

    if (adpcproj_estimate.object$linkfunc == "power5") {
      if (recent) {
        rate <- (agepar +
          driftpar * noobsper +
          driftrecent * driftmp +
          cohpar)^5
      } else {
        rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)^5
      }
    } else if (adpcproj_estimate.object$linkfunc == "log") {
      if (recent) {
        rate <- exp(
          agepar + driftpar * noobsper + driftrecent * driftmp + cohpar
        )
      } else {
        rate <- exp(agepar + driftpar * (noobsper + driftmp) + cohpar)
      }
    } else if (adpcproj_estimate.object$linkfunc == "sqrt") {
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
    noperiod = adpcproj_estimate.object$noperiod,
    gofpvalue = adpcproj_estimate.object$gofpvalue,
    recent = recent,
    pvaluerecent = adpcproj_estimate.object$pvaluerecent,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    distribution = adpcproj_estimate.object$distribution,
    startuseage = startuseage,
    startage = adpcproj_estimate.object$startage,
    glm = adpcproj_estimate.object$glm
  )
  class(res) <- "adpcproj"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' adpcproj_get_predictions
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
#' @keywords internal
adpcproj_get_predictions <- function(
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
    rlang::abort("Variable \"adpcproj.object\" must be of type \"adpcproj\"")
  }

  validate_getpred_inputs(byage, standpop, incidence, agegroups, excludeobs)

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


#' Get ADPCPROJ Projections
#'
#' `get_projections()` extracts annual projection results from an `adpcproj()`
#'   object.
#'
#' @inheritParams canproj
#' @inheritParams adpcproj_get_predictions
#' @param object Output from `adpcproj()`.
#'
#' @return `get_projections()` returns a `data.frame`.
#'
#' @rdname adpcproj
#'
#' @export
get_projections.adpcproj <- function(
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
      "\"cdat\" must match adpcproj.object periods (floor(ncol(cdat) / 5) == noperiod)"
    )
  }

  if (floor(ncol(pdat) / 5) != ncol(object$predictions)) {
    rlang::abort(
      "\"pdat\" must match adpcproj.object periods (floor(ncol(pdat) / 5) == ncol(predictions))"
    )
  }

  r0 <- adpcproj_get_predictions(object, incidence = T)
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
    rlang::abort("Variable \"object\" must be of type \"adpcproj\"")
  }

  validate_summary_inputs(
    printpred = printpred,
    printcall = printcall,
    digits = digits
  )

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
    "First age group estimated (startage):"
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


#' plot.adpcproj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' adpcproj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
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
  ylab = "Rate per 100,000 people",
  main = "",
  lty = c(1, 3),
  col = c("black", "azure4"),
  ...
) {
  if (!inherits(x, "adpcproj")) {
    rlang::abort("Variable \"x\" must be of type \"adpcproj\"")
  }

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

  plot <- ggplot2::ggplot(
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

  return(plot)
}
