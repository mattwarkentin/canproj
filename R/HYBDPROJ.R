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
    nagg <- select_nagg(cdat, pdat, standpop, projfor)
  }

  aggdata <- aggregate_data(cdat, pdat, nagg)
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

  est <- hybdproj_estimate(
    cases = cases,
    pyr = pyr,
    nagg = nagg,
    ncase = ncase,
    linkfunc = linkfunc,
    pD = pD,
    pGOF = pGOF
  )
  pred <- hybdproj_predict(
    hybdproj_estimate.object = est,
    cuttrd = cuttrd,
    shortp = shortp
  )
  return(pred)
}


#' hybdproj_estimate
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
hybdproj_estimate <- function(
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

  model <- select_hybd_method(casesB, apdata, nagg, pGOF, dnoperiods, pD)
  fmodel <- model[[1]]
  apdatan <- model[[2]]

  if (fmodel != "nba-specific") {
    if (fmodel == "common-trend") {
      formula_string <- "Cases ~ as.factor(Age) + Period - 1"
    } else if (fmodel == "age-specific") {
      formula_string <- "Cases ~ as.factor(Age) + as.factor(Age) * Period - Period - 1"
    } else {
      # average model
      formula_string <- "Cases ~ as.factor(Age) - 1"
    }

    res.glm <- get_glm(formula_string, apdatan, linkfunc)
  } else {
    formula_string <- " Cases ~ as.factor(Age) + as.factor(Age) * Period - Period - 1"
    res.glm <- get_glm_nb(formula_string, apdatan, linkfunc)
  }

  pvalue <- 1 - stats::pchisq(res.glm$deviance, res.glm$df.residual)
  res <- list(
    glm = res.glm,
    cases = cases,
    pyr = pyr,
    agrpave = agpave,
    lastper = model[[3]],
    cuty = model[[4]],
    noperiod = model[[3]],
    noyearagg = nagg,
    nocaseagp = ncase,
    linkfunc = linkfunc,
    agrpmod = agpreg,
    projbase = model[[5]],
    finalmod = fmodel,
    gofpvalue = pvalue
  )
  class(res) <- "hybdproj_estimate"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' hybdproj_predict
#'
#' Extrapolate estimated trend from the final model
#'
#' @inheritParams canproj
#'
#' @param hybdproj_estimate.object An object based on the `hybdproj_estimate` function.
#'
#' @return A `list()`.
#'
#' @export
hybdproj_predict <- function(
  hybdproj_estimate.object,
  cuttrd = 0.04,
  shortp = 0
) {
  if (!inherits(hybdproj_estimate.object, "hybdproj_estimate")) {
    stop(
      "Variable \"hybdproj_estimate.object\" must be of type \"hybdproj_estimate\""
    )
  }

  cases <- hybdproj_estimate.object$cases
  pyr <- hybdproj_estimate.object$pyr
  agpreg <- hybdproj_estimate.object$agrpmod
  agpave <- hybdproj_estimate.object$agrpave
  fmodel <- hybdproj_estimate.object$finalmod
  lastper <- hybdproj_estimate.object$lastper

  noobsper <- ncol(cases)
  nototper <- ncol(pyr)
  nonewpred <- nototper - noobsper

  cuttrend <- get_hybd_cuttrend(
    shortp,
    nonewpred,
    hybdproj_estimate.object$noyearagg,
    cuttrd
  )
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

  obsrate <- (as.matrix(cases[, (noobsper - 4):noobsper]) /
    as.matrix(pyr[, (noobsper - 4):noobsper])) |>
    apply(1, mean)

  coef <- hybdproj_estimate.object$glm$coef
  if (fmodel == "age-specific" || fmodel == "nba-specific") {
    m.eff <- cbind(
      coef[1:(length(coef) / 2)],
      coef[(length(coef) / 2 + 1):length(coef)],
      agpreg
    )
  } else if (fmodel == "common-trend") {
    m.eff <- cbind(
      coef[1:(length(coef) - 1)],
      rep(coef[length(coef)], length(agpreg)),
      agpreg
    )
  } else {
    m.eff <- cbind(coef, rep(0, length(agpreg)), agpreg)
  }

  row.names(m.eff) <- NULL
  mcoef <- cbind(rep(NA, length(agpave)), rep(0, length(agpave)), agpave) |>
    rbind(m.eff) |>
    data.frame()
  acoef <- mcoef[with(mcoef, order(mcoef[, 3])), ]
  row.names(acoef) <- NULL
  colnames(acoef) <- c("a.eff", "p.eff", "agrp")

  rate_mat <- matrix(
    acoef$p.eff,
    length(acoef$p.eff),
    length(driftmp)
  ) |>
    sweep(MARGIN = 2, (lastper + driftmp), `*`)

  if (hybdproj_estimate.object$linkfunc == "power5") {
    rate <- (acoef$a.eff + rate_mat)^5
  } else if (hybdproj_estimate.object$linkfunc == "log") {
    rate <- exp(acoef$a.eff + rate_mat)
  } else if (hybdproj_estimate.object$linkfunc == "sqrt") {
    rate <- (acoef$a.eff + rate_mat)^2
  } else {
    rate <- (acoef$a.eff + rate_mat)
  }

  na_groups <- is.na(acoef$a.eff)
  rate[na_groups, ] <- rep(obsrate[na_groups], length(driftmp))

  datatable[, (noobsper + 1):nototper] <- as.data.frame(
    rate *
      pyr[, (noobsper + 1):nototper]
  )

  res <- list(
    predictions = datatable,
    pyr = pyr,
    cuttrd = cuttrd,
    shortp = shortp,
    cuttrend = cuttrend,
    nopred = nonewpred,
    noperiod = hybdproj_estimate.object$noperiod,
    lastperiod = lastper,
    noobsper = noobsper,
    nototper = nototper,
    noyearagg = hybdproj_estimate.object$noyearagg,
    nocaseagp = hybdproj_estimate.object$nocaseagp,
    agrpave = hybdproj_estimate.object$agrpave,
    agrpmod = hybdproj_estimate.object$agrpmod,
    linkfunc = hybdproj_estimate.object$linkfunc,
    projbase = hybdproj_estimate.object$projbase,
    finalmod = hybdproj_estimate.object$finalmod,
    gofpvalue = hybdproj_estimate.object$gofpvalue,
    glm = hybdproj_estimate.object$glm
  )

  class(res) <- "hybdproj"
  attr(res, "Call") <- sys.call()
  return(res)
}


#' hybdproj_get_predictions
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
hybdproj_get_predictions <- function(
  hybdproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
) {
  if (!inherits(hybdproj.object, "hybdproj")) {
    stop("Variable \"hybdproj.object\" must be of type \"hybdproj\"")
  }

  if (missing(byage)) {
    byage <- ifelse(is.null(standpop), T, F)
  }

  validate_getpred_inputs(byage, standpop, incidence, agegroups)

  res <- make_pred_table(
    hybdproj.object,
    agegroups,
    standpop,
    byage,
    incidence,
    excludeobs
  )

  return(res)
}


#' hybdproj_get_projections
#'
#' Extract projection results.
#'
#' @inheritParams hybdproj_get_predictions
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @export
hybdproj_get_projections <- function(
  cdat,
  pdat,
  startp,
  hybdproj.object,
  standpop = NULL,
  Ave5 = FALSE,
  sum5 = TRUE
) {
  if (!is.null(standpop)) {
    S7::check_is_S7(standpop, StandardPopulation)
  }

  finalmod <- hybdproj.object$finalmod
  r0 <- hybdproj_get_predictions(hybdproj.object, incidence = T)
  nagg <- hybdproj.object$noyearagg
  outasp <- interpolate_age_specific_rates(
    r0,
    cdat,
    pdat,
    startp = startp,
    nagg = nagg
  )

  if (finalmod == "average" & Ave5) {
    mod <- ave5proj(cdat, pdat, startp, sum5 = sum5)
    outasp <- mod$agsproj
  }

  if (is.null(standpop)) {
    return(outasp)
  } else {
    outann <- standardize_annual_rates(outasp, pdat, standpop)
    return(outann)
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
#' @inheritParams hybdproj_get_predictions
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

  indat <- hybdproj_get_projections(
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
