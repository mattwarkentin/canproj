#' Data aggregation
#'
#' Aggregation of data by years
#'
#' @inheritParams canproj
#'
#' @return A `list()`.
#'
#' @keywords internal
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
  cases <- as.data.frame(t(stats::aggregate(
    t(cdat0),
    list(Period = pern),
    sum
  )))[-1, ]
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
  pyr1 <- as.data.frame(t(stats::aggregate(
    t(pdat0),
    list(Period = pern),
    sum
  )))[-1, ]

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
  pyr2 <- as.data.frame(t(stats::aggregate(
    t(pdatn),
    list(Period = pernp),
    sum
  )))[-1, ]

  # combine population data:
  pyr <- as.data.frame(cbind(pyr1, pyr2))
  colnames(pyr) <- 1:(nper + nperp)
  rownames(pyr) <- 1:19
  return(list(cases = cases, pyr = pyr))
}
