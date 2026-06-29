#' Convert to age-standardized
#'
#' Calculate age-standardized annual rates and counts
#'
#' @inheritParams canproj
#' @param rates Observed and projected age-specific rates.
#'
#' @keywords internal
standardize_annual_rates <- function(rates, pdat, stdpop) {
  S7::check_is_S7(stdpop, StandardPopulation)
  c1 <- rates * pdat / 100000
  case <- round(colSums(c1), 0)
  case[case < 0] <- 0

  a1 <- rates * stdpop@weights
  asr <- round(colSums(a1), 6)
  asr[asr < 0] <- 0

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
