#' Canada standard population 2011
#'
#' @return A `vector()`.
#'
#' @export
stdpop_Canada_2011 <- function() {
  ca11
}

#' World standard population
#'
#' @return A `vector()`.
#'
#' @export
stdpop_WHO <- function() {
  wdsd
}

# Standard Canadian population, updated June 2019
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
