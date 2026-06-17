#' Observed age-standardized rates
#'
#' Calculate observed age-standardized rates and total numbers
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
obasr <- function(cdat, pdat, stdpop) {
  S7::check_is_S7(stdpop, StandardPopulation)

  numyears <- dim(cdat)[2]
  popu <- pdat[, 1:numyears]
  rr <- matrix(NA, nrow(cdat), numyears)
  for (i in 1:numyears) {
    rr[, i] <- 100000 * cdat[, i] / popu[, i]
  }
  aspr <- rr * stdpop@weights
  asr <- colSums(aspr)
  case <- colSums(cdat)
  return(cbind(asr, case))
}
