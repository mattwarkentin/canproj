#' Average5
#'
#' R functions for projection of cancer incidence/mortality using the average
#'   methods based on the recent 5 years data:
#'   (i) average numbers and population sizes then calculate rates (default), or
#'   (ii) average the calculated yearly rates (sum5=T)
#'
#' @inheritParams canproj
#'
#' @return A `list()`.
#'
#' @export
ave5proj <- function(cdat, pdat, startp, sum5 = TRUE) {
  if (ncol(cdat) < 5) {
    rlang::abort("\"cdat\" must have at least 5 periods")
  }

  if (ncol(pdat) < 10) {
    rlang::abort("\"pdat\" must have at least 10 periods")
  }

  if (!inherits(sum5, "logical")) {
    rlang::abort("\"sum5\" must be of type \"logical\"")
  }

  validate_projection_inputs(
    cdat,
    pdat,
    startp = startp
  )

  nc <- ncol(cdat)
  np <- ncol(pdat)
  noypred <- np - nc
  ngroups <- nrow(cdat)

  obasr <- matrix(NA, ngroups, nc)

  for (i in 1:nc) {
    obasr[, i] <- 100000 * cdat[, i] / pdat[, i]
  }

  datatab <- matrix(NA, ngroups, np)
  datatab[, 1:nc] <- as.matrix(obasr)
  datatab <- data.frame(datatab)

  row.names(datatab) <- 1:ngroups
  colnames(datatab) <- (startp - nc):(startp + (np - nc) - 1)

  if (sum5) {
    obsinc <- rowSums(cdat[, (nc - 4):nc]) / rowSums(pdat[, (nc - 4):nc])
  } else {
    obsinc <- rowSums(cdat[, (nc - 4):nc] / pdat[, (nc - 4):nc]) / 5
  }

  if (sum(is.na(obsinc))) {
    obsinc[is.na(obsinc)] <- 0
  }

  datatab[, (nc + 1):np] <- 100000 * obsinc

  res <- list(
    agsproj = datatab,
    startp = startp,
    sum5 = sum5,
    noypred = noypred
  )

  class(res) <- "ave5proj"

  attr(res, "Call") <- sys.call()

  return(res)
}


#' ave5proj_get_projections
#'
#' Extract the projection results
#'
#' @inheritParams canproj
#' @param ave5proj.object An object based on the `ave5proj()` function.
#'
#' @return A `data.frame()`.
#'
#' @export
ave5proj_get_projections <- function(pdat, ave5proj.object, standpop = NULL) {
  if (!inherits(ave5proj.object, "ave5proj")) {
    rlang::abort("Variable \"ave5proj.object\" must be of type \"ave5proj\"")
  }

  validate_getproj_inputs(
    pdat = pdat,
    standpop = standpop,
    ave5proj.object
  )

  outasp <- ave5proj.object$agsproj

  if (is.null(standpop)) {
    return(outasp)
  } else {
    annproj <- standardize_annual_rates(
      outasp,
      pdat,
      ave5proj.object$startp,
      standpop,
      check = FALSE
    )
    return(annproj)
  }
}


#' summary.ave5proj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the `ave5proj()` function.
#' @param printcall Whether to print function `Call` for ave5proj.object (`T` or `F`).
#' @param ... Other parameters
#'
#' @return An information table describing the method used.
#'
#' @export
summary.ave5proj <- function(object, printcall = FALSE, ...) {
  if (!inherits(object, "ave5proj")) {
    rlang::abort("Variable \"ave5proj.object\" must be of type \"ave5proj\"")
  }

  validate_summary_inputs(printcall = printcall)

  method <- "Five-Year Average"

  if (is.null(object$sum5)) {
    sum5 <- "5-year period"
  } else {
    sum5 <- "average yearly-rates"
  }

  noypred <- object$noypred

  cat("\nPrediction done with:\n")

  moptions <- matrix(NA, 3, 2)

  moptions[, 1] <- c(
    "Method:",
    "Age-Specific Rate by:",
    "Number of prediction years:"
  )

  moptions[, 2] <- c(method, sum5, noypred)

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
