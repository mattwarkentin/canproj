#' Convert to age-standardized
#'
#' Calculate age-standardized annual rates and counts
#'
#' @inheritParams canproj
#' @param rates Observed and projected age-specific rates.
#'
#' @keywords internal
asry <- function(rates, pdat, standpop) {
  c1 <- rates * pdat / 100000
  case <- round(colSums(c1), 0)
  case[case < 0] <- 0

  a1 <- rates * standpop
  asr <- round(colSums(a1), 6)
  asr[asr < 0] <- 0

  return(cbind(asr, case))
}
