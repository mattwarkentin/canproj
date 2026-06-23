#' Get glm
#'
#' Get glm output for given link function
#'
#' @param data Data to use.
#' @param formula_string A string.
#' @inheritParams canproj
#'
#' @keywords internal
get_glm <- function(formula_string, data, linkfunc) {
  if (linkfunc == "power5") {
    y <- data$y
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

    family <- power5link
  } else if (linkfunc == "log") {
    formula_string <- paste0(formula_string, " + offset(log(y))")
    family <- stats::poisson(link = log)
  } else if (linkfunc == "sqrt") {
    formula_string <- sub("Cases", "Cases / y", formula_string)
    family <- stats::poisson(link = sqrt)
  } else if (linkfunc == "identity") {
    formula_string <- sub("Cases", "Cases / y", formula_string)
    family <- stats::poisson(link = identity)
  } else {
    stop("Unknown \"linkfunc\"")
  }

  suppressWarnings(
    res.glm <- stats::glm(
      stats::as.formula(formula_string),
      data = data,
      family = family
    )
  )

  return(res.glm)
}

#' Get glm nb
#'
#' Get glm nb output for given link function
#'
#' @param data Data to use.
#' @param formula_string A string.
#' @inheritParams canproj
#'
#' @keywords internal
get_glm_nb <- function(formula_string, data, linkfunc) {
  if (linkfunc == "power5") {
    y <- data$y
    glmnb <- suppressWarnings(MASS::glm.nb(
      stats::as.formula(paste0(formula_string, " + offset(log(y))")),
      data = data,
      link = log
    ))
    theta <- suppressWarnings(as.numeric(MASS::theta.md(
      data$Cases,
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

    family <- nbpower5link
  } else if (linkfunc == "log") {
    formula_string <- paste0(formula_string, " + offset(log(y))")
    link <- "log"
  } else if (linkfunc == "sqrt") {
    formula_string <- sub("Cases", "Cases / y", formula_string)
    link <- "sqrt"
  } else if (linkfunc == "identity") {
    formula_string <- sub("Cases", "Cases / y", formula_string)
    link <- "identity"
  } else {
    stop("Unknown \"linkfunc\"")
  }

  if (linkfunc == "power5") {
    res.glm <- suppressWarnings(stats::glm(
      stats::as.formula(formula_string),
      data = data,
      family = nbpower5link
    ))
  } else {
    res.glm <- suppressWarnings(do.call(
      MASS::glm.nb,
      list(
        formula = stats::as.formula(formula_string),
        data = data,
        link = link
      )
    ))
  }

  return(res.glm)
}
