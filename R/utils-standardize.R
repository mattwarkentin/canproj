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
  m <- dim(cdat)[2]
  popu <- pdat[, 1:m]
  rr <- matrix(NA, nrow(cdat), m)
  for (i in 1:m) {
    rr[, i] <- 100000 * cdat[, i] / popu[, i]
  }
  aspr <- rr * standpop
  asr <- apply(aspr, 2, sum)
  case <- apply(cdat, 2, sum)
  return(cbind(asr, case))
}
