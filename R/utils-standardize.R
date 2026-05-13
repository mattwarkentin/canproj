#' Observed age-standardized rates
#'
#' Calculate observed age-standardized rates and total numbers
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
obasr <- function(cdat, pdat, standpop = ca11) {
  m <- dim(cdat)[2]
  popu <- pdat[, 1:m]
  rr <- matrix(NA, 19, m)
  for (i in 1:m) {
    rr[, i] <- 100000 * cdat[, i] / popu[, i]
  }
  aspr <- rr * standpop
  asr <- apply(aspr, 2, sum)
  case <- apply(cdat, 2, sum)
  return(cbind(asr, case))
}

#' Age-standardized rates standard error
#'
#' Calculate age-standardized rates and corresponding standard error
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
asrsd <- function(cdat, pdat, standpop = ca11) {
  m <- dim(cdat)[2]
  popu <- pdat[, 1:m]
  rr <- matrix(NA, 19, m)
  ww <- matrix(NA, 19, m)
  for (i in 1:m) {
    rr[, i] <- cdat[, i] / popu[, i]
    ww[, i] <- cdat[, i] * (standpop / popu[, i])^2
  }
  aspr <- rr * standpop
  asr <- 100000 * apply(aspr, 2, sum)
  asd <- 100000 * sqrt(apply(ww, 2, sum))
  case <- apply(cdat, 2, sum)
  return(cbind(asr, asd))
}

#' Percentage change due to risk, population, and aging
#'
#' Calculating percentage changes due to risk, population growth, and aging.
#'
#' @inheritParams canproj
#'
#' @param byear Reference calendar year
#' @param cyear Comparison calendar year
#' @param starty First calendar year in historical data
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
chper <- function(aspr, pdat, byear, cyear, starty = NULL) {
  # Check data:
  if (dim(aspr)[2] != dim(pdat)[2]) {
    stop("\"aspr\" and \"pdat\" must have data for same calendar years")
  }
  # Identify the first calendar year:
  if (is.null(starty)) {
    caly <- as.numeric(names(aspr))
    starty <- caly[1]
  }
  # Identify the numbers of columns for baseline and comparison:
  n1 <- byear - starty + 1
  n2 <- cyear - starty + 1
  asr12 <- aspr[, c(n1, n2)]
  pat12 <- pdat[, c(n1, n2)]
  # Total numbers in baseline and comparison:
  N1 <- sum(pat12[, 1] * asr12[, 1] / 100000)
  N2 <- sum(pat12[, 2] * asr12[, 2] / 100000)
  # Crude rate in baseline year:
  R1 <- N1 / sum(pat12[, 1])
  # Total number in comparison if no change in risk:
  N3 <- sum(pat12[, 2] * asr12[, 1] / 100000)
  # Total number in comparison if no change by aging:
  N4 <- R1 * sum(pat12[, 2])

  # Overall change:
  C.all <- 100 * (N2 - N1) / N1
  # Change due to population growth and aging:
  C.pop <- 100 * (N3 - N1) / N1
  # Change due to risk:
  C.rsk <- C.all - C.pop
  # Change due to population growth:
  C.gwh <- 100 * (N4 - N1) / N1
  # Change due to population aging:
  C.age <- C.pop - C.gwh

  # Output:
  out <- data.frame(matrix(c(N1, N2, C.all, C.rsk, C.gwh, C.age), ncol = 6))
  colnames(out) <- c(
    "ref.case",
    "comp.case",
    "overall",
    "risk",
    "p.growth",
    "p.aging"
  )
  return(out)
}
