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
ave5proj <- function(cdat, pdat, startp, sum5 = NULL) {
  nc <- dim(cdat)[2]
  np <- dim(pdat)[2]
  noypred <- np - nc

  # fill in observed rates:
  obasr <- matrix(NA, 19, nc)

  for (i in 1:nc) {
    obasr[, i] <- 100000 * cdat[, i] / pdat[, i]
  }

  datatab <- matrix(NA, 19, np)

  datatab[, 1:nc] <- as.matrix(obasr)

  datatab <- data.frame(datatab)

  row.names(datatab) <- c(
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

  colnames(datatab) <- (startp - nc):(startp + (np - nc) - 1)

  # Calculate predictions in age-specific rates:
  for (age in 1:19) {
    if (is.null(sum5)) {
      obsinc <- apply(cdat[age, (nc - 5):nc], 1, sum) /
        apply(pdat[age, (nc - 5):nc], 1, sum)
    } else {
      obsinc <- cdat[age, (nc - 5):nc] / pdat[age, (nc - 5):nc]
    }
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0
    }
    if (is.null(sum5)) {
      datatab[age, (nc + 1):np] <- 100000 * obsinc
    } else {
      datatab[age, (nc + 1):np] <- 100000 *
        (obsinc[, 1] + obsinc[, 2] + obsinc[, 3] + obsinc[, 4] + obsinc[, 5]) /
        5
    }
  }

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


#' ave5proj.getproj
#'
#' Extract the projection results
#'
#' @inheritParams canproj
#' @param ave5proj.object An object based on the `ave5proj()` function.
#'
#' @return A `data.frame()`.
#'
#' @export
ave5proj.getproj <- function(pdat, ave5proj.object, standpop = NULL) {
  outasp <- ave5proj.object$agsproj

  if (is.null(standpop)) {
    return(outasp)
  } else {
    annproj <- asry(outasp, pdat, standpop = standpop)
    return(annproj)
  }
}


#' summary.ave5proj
#'
#' Summarize information on projection method used.
#'
#' @param object An object based on the `ave5proj()` function.
#' @param printcall Whether to print function `Call` for ave5proj.object (`T` or `F`).
#' @param digits Number of digits in output, default (`0`) is integer only.
#' @param ... Other parameters
#'
#' @return An information table describing the method used.
#'
#' @export
summary.ave5proj <- function(object, printcall = FALSE, digits = 0, ...) {
  if (!inherits(object, "ave5proj")) {
    stop("Variable \"ave5proj.object\" must be of type \"ave5proj\"")
  }

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
