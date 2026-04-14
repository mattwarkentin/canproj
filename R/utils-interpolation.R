## Linear interpolation:
## Convert to annual age-specific rates
## by linear interpolation for each two-points segment
## if rate is projected age-specific rates in period:
asrpy <- function(rate, cdat, pdat, startp, nagg) {
  # nagg: number of years used for aggregation: 1, 2, ..., 5
  # r0: 19*m matrix, age-specific rate of m periods
  nt <- dim(cdat)[2]
  nperd <- floor(nt / nagg) # of observed periods if aggregated by nagg
  np <- dim(pdat)[2]
  if (nperd <= 25) {
    nypop <- dim(pdat)[2] # of years for observed and projected data
    ny <- dim(cdat)[2] # of years for observed data
    nyp <- nypop - ny # of years for projection
    pdat0 <- pdat
  } else {
    ny <- 25
    nypop <- dim(pdat)[2] - (nt - 25)
    nyp <- nypop - ny
    pdat0 <- pdat[, -c(1:(nt - 25))]
  }

  nper <- floor(ny / nagg) # of observed periods if aggregated by nagg
  nperp <- floor(nyp / nagg) # of projected period if aggregated by nagg

  nycp <- nperp * nagg # of projection years be aggregated
  ry <- nyp - nycp # of rest projection years not aggregated

  pdat1 <- pdat0[, 1:ny] # yearly observed population
  pdat2 <- pdat0[, (ny + 1):nypop] # yearly projected population

  mt <- dim(rate)[2] # of periods from output
  r0 <- as.matrix(rate[, (nper + 1):mt]) # of projected periods
  m <- dim(r0)[2]
  # producing the end periods rates:
  rc1 <- rate[, nper]
  if (m > 1) {
    rc2 <- 2 * r0[, m] - r0[, (m - 1)]
  } else {
    rc2 <- 2 * r0[, m] - rc1
  }
  rc3 <- 2 * rc2 - r0[, m]
  r1 <- cbind(rc1, r0, rc2)
  # producing annual age-specific rates:
  rr <- matrix(NA, 19, nyp)
  if (nagg == 1) {
    rr <- r0
  } else if (nagg == 2) {
    for (i in 2:(m + 1)) {
      rr[, (i - 2) * 2 + 1] <- (1 / 4) * r1[, i - 1] + (3 / 4) * r1[, i]
      rr[, (i - 2) * 2 + 2] <- (1 / 4) * r1[, i + 1] + (3 / 4) * r1[, i]
    }
    if (nyp > nycp) {
      rr[, nyp] <- (1 / 4) * r1[, (m + 1)] + (3 / 4) * r1[, (m + 2)]
    }
  } else if (nagg == 3) {
    for (i in 2:(m + 1)) {
      rr[, (i - 2) * 3 + 1] <- (1 / 3) * r1[, i - 1] + (2 / 3) * r1[, i]
      rr[, (i - 2) * 3 + 2] <- (0 / 3) * r1[, i - 1] + (3 / 3) * r1[, i]
      rr[, (i - 2) * 3 + 3] <- (1 / 3) * r1[, i + 1] + (2 / 3) * r1[, i]
    }
    if (nyp == (nycp + 1)) {
      rr[, nyp] <- (1 / 3) * r1[, (m + 1)] + (2 / 3) * r1[, (m + 2)]
    }
    if (nyp == (nycp + 2)) {
      rr[, (nyp - 1)] <- (1 / 3) * r1[, (m + 1)] + (2 / 3) * r1[, (m + 2)]
      rr[, nyp] <- rc2
    }
  } else if (nagg == 4) {
    for (i in 2:(m + 1)) {
      rr[, (i - 2) * 4 + 1] <- (3 / 8) * r1[, i - 1] + (5 / 8) * r1[, i]
      rr[, (i - 2) * 4 + 2] <- (1 / 8) * r1[, i - 1] + (7 / 8) * r1[, i]
      rr[, (i - 2) * 4 + 3] <- (1 / 8) * r1[, i + 1] + (7 / 8) * r1[, i]
      rr[, (i - 2) * 4 + 4] <- (3 / 8) * r1[, i + 1] + (5 / 8) * r1[, i]
    }
    if (nyp == (nycp + 1)) {
      rr[, nyp] <- (3 / 8) * r1[, (m + 1)] + (5 / 8) * r1[, (m + 2)]
    }
    if (nyp == (nycp + 2)) {
      rr[, (nyp - 1)] <- (3 / 8) * r1[, (m + 1)] + (5 / 8) * r1[, (m + 2)]
      rr[, nyp] <- (1 / 8) * r1[, (m + 1)] + (7 / 8) * r1[, (m + 2)]
    }
    if (nyp == (nycp + 3)) {
      rr[, (nyp - 2)] <- (3 / 8) * r1[, (m + 1)] + (5 / 8) * r1[, (m + 2)]
      rr[, (nyp - 1)] <- (1 / 8) * r1[, (m + 1)] + (7 / 8) * r1[, (m + 2)]
      rr[, nyp] <- (7 / 8) * r1[, (m + 2)] + (1 / 8) * rc3
    }
  } else if (nagg == 5) {
    for (i in 2:(m + 1)) {
      rr[, (i - 2) * 5 + 1] <- (2 / 5) * r1[, i - 1] + (3 / 5) * r1[, i]
      rr[, (i - 2) * 5 + 2] <- (1 / 5) * r1[, i - 1] + (4 / 5) * r1[, i]
      rr[, (i - 2) * 5 + 3] <- (0 / 5) * r1[, i - 1] + (5 / 5) * r1[, i]
      rr[, (i - 2) * 5 + 4] <- (1 / 5) * r1[, i + 1] + (4 / 5) * r1[, i]
      rr[, (i - 2) * 5 + 5] <- (2 / 5) * r1[, i + 1] + (3 / 5) * r1[, i]
    }
    if (nyp == (nycp + 1)) {
      rr[, nyp] <- (2 / 5) * r1[, (m + 1)] + (3 / 5) * r1[, (m + 2)]
    }
    if (nyp == (nycp + 2)) {
      rr[, (nyp - 1)] <- (2 / 5) * r1[, (m + 1)] + (3 / 5) * r1[, (m + 2)]
      rr[, nyp] <- (1 / 5) * r1[, (m + 1)] + (4 / 5) * r1[, (m + 2)]
    }
    if (nyp == (nycp + 3)) {
      rr[, (nyp - 2)] <- (2 / 5) * r1[, (m + 1)] + (3 / 5) * r1[, (m + 2)]
      rr[, (nyp - 1)] <- (1 / 5) * r1[, (m + 1)] + (4 / 5) * r1[, (m + 2)]
      rr[, nyp] <- rc2
    }
    if (nyp == (nycp + 4)) {
      rr[, (nyp - 3)] <- (2 / 5) * r1[, (m + 1)] + (3 / 5) * r1[, (m + 2)]
      rr[, (nyp - 2)] <- (1 / 5) * r1[, (m + 1)] + (4 / 5) * r1[, (m + 2)]
      rr[, (nyp - 1)] <- rc2
      rr[, nyp] <- (1 / 5) * rc3 + (4 / 5) * rc2
    }
  } else {
    stop("Years of aggregation \"nagg\" must be integer between 1 and 5")
  }

  # fill in observed age-specific rates per 100,000:
  obasr <- matrix(NA, 19, nt)
  for (i in 1:nt) {
    obasr[, i] <- 100000 * cdat[, i] / pdat[, i]
  }
  datatab <- matrix(NA, 19, np)
  datatab[, 1:nt] <- as.matrix(obasr)
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
  colnames(datatab) <- (startp - nt):(startp + nyp - 1)
  # fill in projected age-specific rates per 100,000:
  datatab[, (nt + 1):np] <- rr
  return(datatab)
}
