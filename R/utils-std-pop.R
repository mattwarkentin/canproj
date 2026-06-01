#' 2011 Canadian standard population, updated June 2019
#'
#' @return A `vector()`.
#'
#' @export
stdpop_Canada_2011 <- function() {
  c(
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
}

#' WHO world standard population, updated June 2026
#'
#' @return A `vector()`.
#'
#' @export
stdpop_WHO <- function() {
  c(
    0.08857,
    0.08687,
    0.08597,
    0.08467,
    0.08217,
    0.07927,
    0.07607,
    0.07148,
    0.06588,
    0.06038,
    0.05368,
    0.04548,
    0.03719,
    0.02959,
    0.02209,
    0.01519,
    0.00910,
    0.00440,
    0.00195
  )
  # Age group 35-39 was increased by 0.00001 so the total sum is 1
  # This was the closes number to rounding up that did not
}
