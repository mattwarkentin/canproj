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
  num_obs_year <- ncol(cdat)
  nagg_obs_period <- floor(num_obs_year / nagg)

  if (num_obs_year == nagg_obs_period * nagg) {
    trunc_cdat <- cdat
  } else {
    trunc_cdat <- cdat[, -c(1:(num_obs_year - nagg_obs_period * nagg))]
  }

  per_groups <- sort(rep(1:nagg_obs_period, nagg))

  cases <- as.data.frame(t(rowsum(t(trunc_cdat), per_groups)))
  colnames(cases) <- 1:nagg_obs_period
  rownames(cases) <- 1:nrow(cdat)

  obs_pop <- pdat[, 1:num_obs_year]
  if (num_obs_year == nagg_obs_period * nagg) {
    trunc_obs_pop <- obs_pop
  } else {
    trunc_obs_pop <- obs_pop[, -c(1:(num_obs_year - nagg_obs_period * nagg))]
  }

  obs_pyr <- as.data.frame(t(rowsum(t(trunc_obs_pop), per_groups)))

  proj_pop <- pdat[, c((num_obs_year + 1):ncol(pdat))]
  num_proj_year <- ncol(pdat) - num_obs_year
  nagg_proj_period <- floor(num_proj_year / nagg)

  if (num_proj_year == nagg_proj_period * nagg) {
    trunc_proj_pop <- proj_pop
  } else {
    trunc_proj_pop <- proj_pop[,
      -c((nagg_proj_period * nagg + 1):num_proj_year)
    ]
  }

  proj_groups <- sort(rep(1:nagg_proj_period, nagg))
  proj_pyr <- as.data.frame(t(rowsum(t(trunc_proj_pop), proj_groups)))

  pyr <- as.data.frame(cbind(obs_pyr, proj_pyr))
  colnames(pyr) <- 1:(nagg_obs_period + nagg_proj_period)
  rownames(pyr) <- 1:nrow(cdat)

  return(list(cases = cases, pyr = pyr))
}
