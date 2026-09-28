#' Project Using Constant Rates
#'
#' Project future cancer cases by applying age-specific cancer rates from the
#'   most recent year of observed data to projected future population
#'   distributions. The function assumes that the most recent observed cancer
#'   rates remain constant over the projection period and applies these rates
#'   to the corresponding age strata in each future year.
#'
#' @inheritParams canproj
#'
#' @return A `list()`.
#'
#' @export
project_constant <- function(cdat, pdat, startp) {
  validate_project_constant_inputs(cdat, pdat, startp)

  cdat <- as.matrix(cdat)
  pdat <- as.matrix(pdat)

  startyear <- startp - ncol(cdat)

  pdat_years <- seq(0, ncol(pdat) - 1, 1) + startyear

  lastknown_cases <- cdat[, ncol(cdat)]

  lastknown_pop <- pdat[, ncol(cdat)]

  constant_aspcr <- lastknown_cases / lastknown_pop

  proj_cases <- round(constant_aspcr * pdat[, (ncol(cdat) + 1):ncol(pdat)])

  datatab <- cbind(cdat, proj_cases)

  colnames(datatab) <- pdat_years

  datatab <- as.data.frame(datatab)

  noypred <- ncol(pdat) - ncol(cdat)

  res <- list(
    agsproj = datatab,
    noypred = noypred,
    cdat = cdat,
    pdat = pdat,
    startp = startp
  )

  class(res) <- c("constant", "proj_model")
  attr(res, "Call") <- sys.call()

  res
}
