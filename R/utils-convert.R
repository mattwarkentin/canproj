#' Convert to age-standardized
#'
#' Calculate age-standardized annual rates and counts
#'
#' @inheritParams canproj
#' @param rr Observed and projected age-specific rates.
#'
#' @keywords internal
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

# asrny <- function(r0, cdat, pdat, startp, nagg, standpop = ca91) {
#   rr <- asrpy(r0, cdat, pdat, startp = startp, nagg)
#   c1 <- rr * pdat / 100000
#   a1 <- rr * standpop
#   asr <- round(apply(a1, 2, sum), 6)
#   case <- round(apply(c1, 2, sum), 0)
#   for (i in 1:length(asr)) {
#     if (asr[i] < 0) asr[i] <- 0
#   }
#   for (i in 1:length(case)) {
#     if (case[i] < 0) case[i] <- 0
#   }
#   return(cbind(asr, case))
# }
