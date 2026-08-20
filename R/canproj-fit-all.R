#' Fit all Canproj models
#'
#' Fits all of the available `canproj()` methods. This function accepts all the
#'   same arguments as `canproj()` and they are passed on to each function call.
#'   The methods fit by this model include: `"nordpred"`, `"adpc-nb"`,
#'   `"ac-poi"`, `"ac-nb"`, `"age-trd-poi"`, `"age-trd-nb"`, `"com-trd"`,
#'   `"age-only"`, and `"ave5"`.
#'
#' @inheritParams canproj
#'
#' @return A `list`. The `list` contains the following items:
#'   - `<selected_method>`: Projection method selected by `canproj()`.
#'   - `nordpred`: `canproj` object produced by the `"nordpred"` method.
#'   - `adpc-nb`: `canproj` object produced by the `"adpc-nb"` method.
#'   - `ac-poi`: `canproj` object produced by the `"ac-poi"` method.
#'   - `ac-nb`:`canproj` object produced by the `"ac-nb"` method.
#'   - `a-s-nb`: `canproj` object produced by the `"a-s-nb"` method.
#'   - `a-s-poi`: `canproj` object produced by the `"a-s-poi"` method.
#'   - `c-t`: `canproj` object produced by the `"c-t"` method.
#'   - `average`: `canproj` object produced by the `"average"` method.
#'   - `ave5`: `canproj` object produced by the `"ave5"` method.
#'
#' @export
canproj_all_methods <- function(
  cdat,
  pdat,
  startp,
  standpop,
  projfor = "incidence",
  nagg = NULL,
  ncase = NULL,
  startage = NULL,
  newcohort = FALSE,
  ave5 = FALSE,
  sum5 = TRUE,
  linkfunc = "power5",
  cuttrd = 0.04,
  shortp = 0,
  pD = 0.05,
  pGOF = 0.05
) {
  methods <- c(
    "nordpred",
    "adpc-nb",
    "ac-poi",
    "ac-nb",
    "age-trd-poi",
    "age-trd-nb",
    "com-trd",
    "age-only",
    "ave5"
  )
  out <- list()

  selected <- canproj(
    cdat,
    pdat,
    startp,
    standpop,
    projfor = projfor,
    nagg = nagg,
    ncase = ncase,
    startage = startage,
    newcohort = newcohort,
    ave5 = ave5,
    sum5 = sum5,
    linkfunc = linkfunc,
    cuttrd = cuttrd,
    shortp = shortp,
    pD = pD,
    pGOF = pGOF
  )

  if (selected$method == "ADPC") {
    if (selected$out$distribution == "Poisson") {
      out$selected_method <- "nordpred"
    } else {
      out$selected_method <- "adpc_nb"
    }
  } else if (selected$method == "AC") {
    if (selected$out$distribution == "Poisson") {
      out$selected_method <- "ac_poi"
    } else {
      out$selected_method <- "ac_nb"
    }
  } else if (selected$method == "Hybrid") {
    if (selected$out$finalmod == "age-specific") {
      out$selected_method <- "age_trd_poi"
    } else if (selected$out$finalmod == "nba-specific") {
      out$selected_method <- "age_trd_nb"
    } else if (selected$out$finalmod == "common-trend") {
      out$selected_method <- "com_trd"
    } else {
      out$selected_method <- "age_only"
    }
  } else {
    # ave5
    out$selected_method <- "ave5"
  }

  for (m in methods) {
    out[[gsub("-", "_", m)]] <- canproj(
      cdat,
      pdat,
      startp,
      standpop,
      projfor = projfor,
      nagg = nagg,
      ncase = ncase,
      startage = startage,
      newcohort = newcohort,
      ave5 = ave5,
      sum5 = sum5,
      methods = m,
      linkfunc = linkfunc,
      cuttrd = cuttrd,
      shortp = shortp,
      pD = pD,
      pGOF = pGOF
    )
  }

  return(out)
}
