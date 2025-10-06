# Data aggregation function:
##aggregating cases and population by nagg=x years
datagg <- function(cdat, pdat, nagg) {
  # aggregating cancer data:
  ny <- dim(cdat)[2] # of years for observed data
  nper <- floor(ny / nagg) # of period if aggregated by nagg
  if (ny == nper * nagg) {
    cdat0 <- cdat
  } else {
    # truncate the beginning years less than nagg
    cdat0 <- cdat[, -c(1:(ny - nper * nagg))]
  }
  pern <- sort(rep(1:nper, nagg)) # period vector
  cases <- as.data.frame(t(aggregate(t(cdat0), list(Period = pern), sum)))[-1, ]
  colnames(cases) <- 1:nper
  rownames(cases) <- 1:19

  # aggregating population data:
  pdat1 <- pdat[, 1:ny] # observed population
  if (ny == nper * nagg) {
    pdat0 <- pdat1
  } else {
    # truncate the beginning years less than nagg
    pdat0 <- pdat1[, -c(1:(ny - nper * nagg))]
  }
  pyr1 <- as.data.frame(t(aggregate(t(pdat0), list(Period = pern), sum)))[-1, ]

  # aggregating projection population data:
  nypop <- dim(pdat)[2]
  pdat2 <- pdat[, c((ny + 1):nypop)] # projected population
  nyp <- nypop - ny
  nperp <- floor(nyp / nagg) # of period if aggregated by nagg
  if (nyp == nperp * nagg) {
    pdatn <- pdat2
  } else {
    # truncate the end years great than nagg
    pdatn <- pdat2[, -c((nperp * nagg + 1):nyp)]
  }
  pernp <- sort(rep(1:nperp, nagg)) # period vector
  pyr2 <- as.data.frame(t(aggregate(t(pdatn), list(Period = pernp), sum)))[-1, ]

  # combine population data:
  pyr <- as.data.frame(cbind(pyr1, pyr2))
  colnames(pyr) <- 1:(nper + nperp)
  rownames(pyr) <- 1:19
  return(list(cases = cases, pyr = pyr))
}


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


## Convert to annual ASRs and Counts:
#asrny <- function(r0, cdat, pdat, startp, nagg, standpop=ca91) {
#    rr <- asrpy(r0, cdat, pdat, startp=startp, nagg)
#    c1 <- rr * pdat / 100000
#    a1 <- rr * standpop
#    asr <- round(apply(a1, 2, sum), 6)
#    case <- round(apply(c1, 2, sum), 0)
#    for (i in 1: length(asr)) {if (asr[i] < 0) asr[i] <- 0}
#    for (i in 1: length(case)) {if (case[i] < 0) case[i] <- 0}
#    return(cbind(asr, case))
#}
## Or
asry <- function(rr, pdat, standpop = ca11) {
  c1 <- rr * pdat / 100000
  a1 <- rr * standpop
  asr <- round(apply(a1, 2, sum), 6)
  case <- round(apply(c1, 2, sum), 0)
  for (i in 1:length(asr)) {
    if (asr[i] < 0) asr[i] <- 0
  }
  for (i in 1:length(case)) {
    if (case[i] < 0) case[i] <- 0
  }
  return(cbind(asr, case))
}


## Standardizing observed rates and total numbers:
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

asrsd <- function(cdat, pdat, standpop = ca11) {
  ## calculate the age-standardized rate with correspond standard error
  ## cdat, pdat: 19 age groups by N years
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


chper <- function(aspr, pdat, byear, cyear, starty = NULL) {
  ## Calculating the percentage changes due to risk, population growth and aging.
  ## pdat: population data over historical and projection years,
  ## aspr: data frame for annual age-specific rates,
  ## byear: baseline calendar year, cyear: comparison calendar year,
  ## starty: the first calendar year in historical data.

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


# Projection plot contains both ASRs and total numbers:
projplot <- function(
  site.asr,
  sex = NULL,
  ma = 2,
  mr = 1.02,
  starty = 1986,
  startp = 2011
) {
  case <- site.asr[, 2]
  asr <- site.asr[, 1]
  n <- dim(site.asr)[1] #total number of years
  n1 <- startp - starty #n1 in observations
  n2 <- n - n1 #n2 in projections
  maxn <- max(case)
  maxr <- max(asr)
  if (is.null(sex)) {
    mycol <- c("red4", "red2", "red4", "red2")
  } else if (sex == "M") {
    mycol <- c("darkblue", "blue", "darkblue", "blue2")
  } else {
    mycol <- c("hotpink", "lightpink", "hotpink", "lightpink2")
  }
  par(mar = c(4, 4, 2, 4), mgp = c(3, 1, 0), cex = 0.7)
  barplot(
    c(case[1:n1], rep(0, n2)),
    col = mycol[1],
    space = 1,
    ylim = c(0, ma * maxn),
    #            cex.axis=1,xlab="", las=1, xaxs="r", border = NA, axis.lty = 1)
    cex.axis = 1,
    xlab = "",
    las = 1,
    xaxs = "r",
    border = NA
  )
  abline(h = 0)
  par(mgp = c(3, 1, 0))
  barplot(
    c(rep(0, n1), case[(n1 + 1):n]),
    add = T,
    col = mycol[2],
    space = 1,
    cex.axis = 1,
    las = 1,
    xaxs = "r",
    border = NA
  )
  par(new = T, mar = c(4, 4, 2, 4), mgp = c(3, 1, 0), cex = 1)
  plot(
    c(starty:(starty + n - 1)),
    c(asr[1:n1], rep("", n2)),
    type = "o",
    pch = 20,
    ylim = c(0, maxr),
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    lwd = 3,
    ylab = "",
    col = mycol[3],
    xlab = "",
    cex.axis = 1
  )
  lines(
    c(starty:(starty + n - 1)),
    c(rep("", (n1 - 1)), asr[n1:n]),
    type = "l",
    pch = 20,
    ylim = c(0, maxr),
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    lwd = 3,
    lty = 2,
    col = mycol[4],
    cex.axis = 1
  )
  axis(4, at = seq(0, (mr * maxr), ceiling(maxr / 6)), las = 1, cex.axis = 1)
  #    axis(4, las=1, cex.axis=1)
  mtext(side = 4, "", line = 2.5, cex = 1)
}


## Standard populations:
# ca91 <- c(0.069465, 0.069454, 0.068034, 0.068495, 0.075016, 0.089944, 0.092400, 0.083388,
#           0.076063, 0.059536, 0.047649, 0.044041, 0.042326, 0.038570, 0.029660, 0.022127,
#           0.013595, 0.010237)
#
# ca96 <- c(0.066235, 0.067985, 0.067716, 0.067841, 0.067761, 0.072914, 0.087030, 0.088510,
#           0.080055, 0.071847, 0.055812, 0.044869, 0.040705, 0.037858, 0.032589, 0.023232,
#           0.015424, 0.011617)

#ca11 <- c(0.055297, 0.052717, 0.055853, 0.065194, 0.068555, 0.069006, 0.067786, 0.066188,
#           0.069474, 0.079199, 0.078365, 0.068518, 0.059705, 0.044636, 0.033597, 0.026769,
#           0.020416, 0.018725)
#Updated June 2019
ca11 <- c(
  0.055325,
  0.052736,
  0.055857,
  0.065121,
  0.068506,
  0.068967,
  0.067743,
  0.066176,
  0.069477,
  0.079209,
  0.078366,
  0.068519,
  0.059696,
  0.044613,
  0.033573,
  0.026761,
  0.020406,
  0.012420,
  0.006529
)

# Segi world population
# wdsd <- c(.12, .1, .09, .09, .08, .08, .06, .06, .06, .06, .05, .04, .04, .03, .02, .01,
#           .005, .005)

# New WHO world population
wdsd <- c(
  0.08855,
  0.08687,
  0.08597,
  0.08474,
  0.08222,
  0.07928,
  0.07605,
  0.07145,
  0.0659,
  0.06038,
  0.05371,
  0.04547,
  0.03723,
  0.02955,
  0.0221,
  0.01515,
  0.00905,
  0.0044,
  0.00193
)
