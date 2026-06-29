#' Linear interpolation
#'
#' Convert to annual age-specific rates by linear interpolation
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
interpolate_age_specific_rates <- function(rate, cdat, pdat, startp, nagg) {
  ngroups <- nrow(cdat)
  tot_obs_years <- ncol(cdat)
  tot_years <- ncol(pdat)

  if (floor(tot_obs_years / nagg) <= 25) {
    ny_pop <- ncol(pdat)
    ny_obs <- ncol(cdat)
    ny_proj <- ny_pop - ny_obs
    pdat_effective <- pdat
  } else {
    ny_obs <- 25
    ny_pop <- ncol(pdat) - (tot_obs_years - 25)
    ny_proj <- ny_pop - ny_obs
    pdat_effective <- pdat[, -c(1:(tot_obs_years - 25))]
  }

  nagg_obs_per <- floor(ny_obs / nagg)
  nagg_proj_year <- floor(ny_proj / nagg) * nagg

  pred_rates <- as.matrix(rate[, (nagg_obs_per + 1):ncol(rate)])
  npred <- dim(pred_rates)[2]

  last_obs_rate <- rate[, nagg_obs_per]
  if (npred > 1) {
    next_rate <- 2 * pred_rates[, npred] - pred_rates[, (npred - 1)]
  } else {
    next_rate <- 2 * pred_rates[, npred] - last_obs_rate
  }
  last_rate <- 2 * next_rate - pred_rates[, npred]
  r1 <- cbind(last_obs_rate, pred_rates, next_rate)

  annual_rates <- interpolate(
    nagg,
    matrix(NA, ngroups, ny_proj),
    r1,
    next_rate,
    last_rate,
    npred,
    nagg_proj_year,
    ny_proj
  )

  obasr <- matrix(NA, ngroups, tot_obs_years)

  for (i in 1:tot_obs_years) {
    obasr[, i] <- 100000 * cdat[, i] / pdat[, i]
  }

  datatab <- matrix(NA, ngroups, tot_years)
  datatab[, 1:tot_obs_years] <- as.matrix(obasr)
  datatab <- data.frame(datatab)
  row.names(datatab) <- 1:ngroups
  colnames(datatab) <- (startp - tot_obs_years):(startp + ny_proj - 1)

  datatab[, (tot_obs_years + 1):tot_years] <- annual_rates
  return(datatab)
}

#' Interpolation
#'
#' Calculate rates by linear interpolation
#'
#' @inheritParams canproj
#' @param matrix Base matrix to work off of.
#' @param rates Matrix of rates for predicted periods.
#' @param next_rate Calculated rate for extra years.
#' @param last_rate Calculated rate for extra years.
#' @param npred Number of predicted periods after aggregation.
#' @param nagg_proj_year Number of aggregated projection years.
#' @param ny_proj Number of years for preojection, not aggregated.
#'
#' @return A `matrix()`.
#'
#' @keywords internal
interpolate <- function(
  nagg,
  matrix,
  rates,
  next_rate,
  last_rate,
  npred,
  nagg_proj_year,
  ny_proj
) {
  if (nagg == 1) {
    matrix <- rates[, -c(1, ncol(rates))]
  } else if (nagg == 2) {
    for (i in 2:(npred + 1)) {
      matrix[, (i - 2) * 2 + 1] <- (1 / 4) *
        rates[, i - 1] +
        (3 / 4) * rates[, i]
      matrix[, (i - 2) * 2 + 2] <- (1 / 4) *
        rates[, i + 1] +
        (3 / 4) * rates[, i]
    }
    if (ny_proj > nagg_proj_year) {
      matrix[, ny_proj] <- (1 / 4) *
        rates[, (npred + 1)] +
        (3 / 4) * rates[, (npred + 2)]
    }
  } else if (nagg == 3) {
    for (i in 2:(npred + 1)) {
      matrix[, (i - 2) * 3 + 1] <- (1 / 3) *
        rates[, i - 1] +
        (2 / 3) * rates[, i]
      matrix[, (i - 2) * 3 + 2] <- (0 / 3) *
        rates[, i - 1] +
        (3 / 3) * rates[, i]
      matrix[, (i - 2) * 3 + 3] <- (1 / 3) *
        rates[, i + 1] +
        (2 / 3) * rates[, i]
    }
    if (ny_proj == (nagg_proj_year + 1)) {
      matrix[, ny_proj] <- (1 / 3) *
        rates[, (npred + 1)] +
        (2 / 3) * rates[, (npred + 2)]
    }
    if (ny_proj == (nagg_proj_year + 2)) {
      matrix[, (ny_proj - 1)] <- (1 / 3) *
        rates[, (npred + 1)] +
        (2 / 3) * rates[, (npred + 2)]
      matrix[, ny_proj] <- next_rate
    }
  } else if (nagg == 4) {
    for (i in 2:(npred + 1)) {
      matrix[, (i - 2) * 4 + 1] <- (3 / 8) *
        rates[, i - 1] +
        (5 / 8) * rates[, i]
      matrix[, (i - 2) * 4 + 2] <- (1 / 8) *
        rates[, i - 1] +
        (7 / 8) * rates[, i]
      matrix[, (i - 2) * 4 + 3] <- (1 / 8) *
        rates[, i + 1] +
        (7 / 8) * rates[, i]
      matrix[, (i - 2) * 4 + 4] <- (3 / 8) *
        rates[, i + 1] +
        (5 / 8) * rates[, i]
    }
    if (ny_proj == (nagg_proj_year + 1)) {
      matrix[, ny_proj] <- (3 / 8) *
        rates[, (npred + 1)] +
        (5 / 8) * rates[, (npred + 2)]
    }
    if (ny_proj == (nagg_proj_year + 2)) {
      matrix[, (ny_proj - 1)] <- (3 / 8) *
        rates[, (npred + 1)] +
        (5 / 8) * rates[, (npred + 2)]
      matrix[, ny_proj] <- (1 / 8) *
        rates[, (npred + 1)] +
        (7 / 8) * rates[, (npred + 2)]
    }
    if (ny_proj == (nagg_proj_year + 3)) {
      matrix[, (ny_proj - 2)] <- (3 / 8) *
        rates[, (npred + 1)] +
        (5 / 8) * rates[, (npred + 2)]
      matrix[, (ny_proj - 1)] <- (1 / 8) *
        rates[, (npred + 1)] +
        (7 / 8) * rates[, (npred + 2)]
      matrix[, ny_proj] <- (7 / 8) * rates[, (npred + 2)] + (1 / 8) * next_rate
    }
  } else if (nagg == 5) {
    for (i in 2:(npred + 1)) {
      matrix[, (i - 2) * 5 + 1] <- (2 / 5) *
        rates[, i - 1] +
        (3 / 5) * rates[, i]
      matrix[, (i - 2) * 5 + 2] <- (1 / 5) *
        rates[, i - 1] +
        (4 / 5) * rates[, i]
      matrix[, (i - 2) * 5 + 3] <- (0 / 5) *
        rates[, i - 1] +
        (5 / 5) * rates[, i]
      matrix[, (i - 2) * 5 + 4] <- (1 / 5) *
        rates[, i + 1] +
        (4 / 5) * rates[, i]
      matrix[, (i - 2) * 5 + 5] <- (2 / 5) *
        rates[, i + 1] +
        (3 / 5) * rates[, i]
    }
    if (ny_proj == (nagg_proj_year + 1)) {
      matrix[, ny_proj] <- (2 / 5) *
        rates[, (npred + 1)] +
        (3 / 5) * rates[, (npred + 2)]
    }
    if (ny_proj == (nagg_proj_year + 2)) {
      matrix[, (ny_proj - 1)] <- (2 / 5) *
        rates[, (npred + 1)] +
        (3 / 5) * rates[, (npred + 2)]
      matrix[, ny_proj] <- (1 / 5) *
        rates[, (npred + 1)] +
        (4 / 5) * rates[, (npred + 2)]
    }
    if (ny_proj == (nagg_proj_year + 3)) {
      matrix[, (ny_proj - 2)] <- (2 / 5) *
        rates[, (npred + 1)] +
        (3 / 5) * rates[, (npred + 2)]
      matrix[, (ny_proj - 1)] <- (1 / 5) *
        rates[, (npred + 1)] +
        (4 / 5) * rates[, (npred + 2)]
      matrix[, ny_proj] <- next_rate
    }
    if (ny_proj == (nagg_proj_year + 4)) {
      matrix[, (ny_proj - 3)] <- (2 / 5) *
        rates[, (npred + 1)] +
        (3 / 5) * rates[, (npred + 2)]
      matrix[, (ny_proj - 2)] <- (1 / 5) *
        rates[, (npred + 1)] +
        (4 / 5) * rates[, (npred + 2)]
      matrix[, (ny_proj - 1)] <- next_rate
      matrix[, ny_proj] <- (1 / 5) * next_rate + (4 / 5) * next_rate
    }
  } else {
    stop("Years of aggregation \"nagg\" must be integer between 1 and 5")
  }

  return(matrix)
}
