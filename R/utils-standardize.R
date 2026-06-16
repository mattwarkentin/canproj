#' Observed age-standardized rates
#'
#' Calculate observed age-standardized rates and total numbers
#'
#' @inheritParams canproj
#'
#' @return A `data.frame()`.
#'
#' @keywords internal
obasr <- function(cdat, pdat, standpop) {
  numyears <- dim(cdat)[2]
  popu <- pdat[, 1:numyears]
  rr <- matrix(NA, nrow(cdat), numyears)
  for (i in 1:numyears) {
    rr[, i] <- 100000 * cdat[, i] / popu[, i]
  }
  aspr <- rr * standpop
  asr <- colSums(aspr)
  case <- colSums(cdat)
  return(cbind(asr, case))
}
