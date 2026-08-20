#' Average5: Five-Year Average Projections
#'
#' R functions for projection of cancer incidence/mortality using the average
#'   methods based on the recent 5 years data:
#'   (i) average numbers and population sizes then calculate rates (default), or
#'   (ii) average the calculated yearly rates (sum5=T)
#'
#' @inheritParams canproj
#'
#' @return `ave5proj()` returns a `list`.
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


#' Get Average5 Projections
#'
#' `get_projections()` extracts the projection results from an `ave5proj()`
#'   object.
#'
#' @inheritParams canproj
#' @param object Output object from `ave5proj()`.
#'
#' @return `get_projections()` returns a `data.frame`.
#'
#' @rdname ave5proj
#' @export
get_projections.ave5proj <- function(object, ..., pdat, standpop = NULL) {
  rlang::check_dots_empty()
  validate_getproj_inputs(
    object = object,
    pdat = pdat,
    standpop = standpop
  )

  outasp <- object$agsproj

  if (is.null(standpop)) {
    return(outasp)
  } else {
    annproj <- standardize_annual_rates(
      outasp,
      pdat,
      object$startp,
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

#' plot.ave5proj
#'
#' Create the graph of the observed and projected age-standardized rates from a
#' ave5proj object
#'
#' @inheritParams plot.canproj
#' @inheritParams canproj
#'
#' @export
plot.ave5proj <- function(
  x,
  pdat,
  standpop,
  startplot = 1,
  xlab = "Calendar Year",
  ylab = "Rate per 100,000 people",
  main = "",
  lty = c(1, 3),
  col = c("black", "azure4"),
  ...
) {
  if (!inherits(x, "ave5proj")) {
    rlang::abort("Variable \"x\" must be of type \"ave5proj\"")
  }
  S7::check_is_S7(standpop, StandardPopulation)

  indat <- get_projections(object = x, pdat = pdat, standpop = standpop)
  indata <- indat[, 1]

  data <- as.data.frame(indata)
  data$year <- as.numeric(names(indata))

  data$Period <- "Observed"
  data$Period[(length(indata) - x$noypred + 1):length(indata)] <- "Projected"
  dupe <- data[(length(indata) - x$noypred), 1:2]
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
