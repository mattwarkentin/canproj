#' Observed age-standardized rates
#'
#' Calculate observed age-standardized rates and total numbers
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
obasr <- function(cdat, pdat, stdpop) {
  S7::check_is_S7(stdpop, StandardPopulation)

  numyears <- dim(cdat)[2]
  popu <- pdat[, 1:numyears]
  rr <- matrix(NA, nrow(cdat), numyears)
  for (i in 1:numyears) {
    rr[, i] <- 100000 * cdat[, i] / popu[, i]
  }
  aspr <- rr * stdpop@weights
  asr <- colSums(aspr)
  case <- colSums(cdat)
  return(cbind(asr, case))
}

#' Age-standardized rates and standard error
#'
#' Calculate age-standardized rates and standard error
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
standardize_rates <- function(cdat, pdat, stdpop = stdpop_Canada_2021) {
  num_years <- ncol(cdat)
  num_agegps <- nrow(cdat)

  pop_data <- pdat[, 1:num_years]

  rr <- matrix(NA, num_agegps, num_years)

  ww <- matrix(NA, num_agegps, num_years)

  for (year in 1:num_years) {
    rr[, year] <- cdat[, year] / pop_data[, year]
    ww[, year] <- cdat[, year] * (stdpop / pop_data[, year])^2
  }

  aspr <- rr * stdpop

  asr <- 100000 * apply(aspr, 2, sum)
  asd <- 100000 * sqrt(apply(ww, 2, sum))

  case <- apply(cdat, 2, sum)

  cbind(asr, asd, case)
}
