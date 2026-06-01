standardize_rates <- function(cdat, pdat, stdpop = stdpop_Canada_2021()) {
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
