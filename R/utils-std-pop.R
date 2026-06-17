#' Standard Population
#'
#' An `S7` class for creating standard populations to be used in other functions
#'   in this package. For convenience, we provide three commonly used standard
#'   populations for Canada (2011 and 2021) and the world (WHO 2000-2025). If
#'   a user requires a custom standard population, you must use
#'   `standard_population()` to construct the standard population and pass the
#'   object as an argument to other functions, as needed.
#'
#' @param name Descriptive name of the standard population.
#' @param strata Character vector of strata labels.
#' @param weights Numeric vector of strata weights.
#' @param metadata Optional. We provide this `metadata` property so that users
#'   can include arbitrary metadata in an object.
#'
#' @return An `S7` `standard_population` object.
#'
#' @import S7
#'
#' @export
standard_population <-
  S7::new_class(
    name = "standard_population",
    package = "canproj",
    properties = list(
      name = S7::new_property(
        class = S7::class_character,
        default = quote(rlang::abort("@name is requied."))
      ),
      strata = S7::new_property(
        class = S7::class_character,
        default = quote(rlang::abort("@strata is requied."))
      ),
      weights = S7::new_property(
        class = S7::class_double,
        default = quote(rlang::abort("@weights is requied.")),
        validator = function(value) {
          if (rlang::is_true(all.equal(sum(value), 1, tolerance = 0.01))) {
            return(NULL)
          }
          "Must sum to 1."
        }
      ),
      metadata = S7::class_list
    ),
    validator = function(self) {
      # Validate `name` ----
      if (length(self@name) > 1) {
        "`name` must have length 1."
        # Validate `weights` and `strata` ----
      } else if (length(self@weights) != length(self@strata)) {
        "`strata` and `weights` must be the same length."
      }
    }
  )

#' @rdname standard_population
#' @export
stdpop_Canada_2021 <-
  standard_population(
    name = "Canadian Standard Population (2021)",
    strata = c(
      "0 to 4 years",
      "5 to 9 years",
      "10 to 14 years",
      "15 to 19 years",
      "20 to 24 years",
      "25 to 29 years",
      "30 to 34 years",
      "35 to 39 years",
      "40 to 44 years",
      "45 to 49 years",
      "50 to 54 years",
      "55 to 59 years",
      "60 to 64 years",
      "65 to 69 years",
      "70 to 74 years",
      "75 to 79 years",
      "80 to 84 years",
      "85 to 89 years",
      "90 years and over"
    ),
    weights = c(
      0.049762,
      0.054199,
      0.055283,
      0.053870,
      0.062877,
      0.069944,
      0.070712,
      0.069120,
      0.065514,
      0.062285,
      0.063488,
      0.070359,
      0.068330,
      0.058224,
      0.048268,
      0.033273,
      0.022004,
      0.013635,
      0.008853
    ),
    metadata = list(
      last_updated = "June 2026"
    )
  )

#' @rdname standard_population
#' @export
stdpop_Canada_2011 <-
  standard_population(
    name = "Canadian Standard Population (2011)",
    strata = c(
      "0 to 4 years",
      "5 to 9 years",
      "10 to 14 years",
      "15 to 19 years",
      "20 to 24 years",
      "25 to 29 years",
      "30 to 34 years",
      "35 to 39 years",
      "40 to 44 years",
      "45 to 49 years",
      "50 to 54 years",
      "55 to 59 years",
      "60 to 64 years",
      "65 to 69 years",
      "70 to 74 years",
      "75 to 79 years",
      "80 to 84 years",
      "85 to 89 years",
      "90 years and over"
    ),
    weights = c(
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
    ),
    metadata = list(
      last_updated = "June 2019"
    )
  )

#' @rdname standard_population
#' @export
stdpop_WHO_2000_2025 <-
  standard_population(
    name = "World (WHO 2000-2025) Standard Population",
    strata = c(
      "0 to 4 years",
      "5 to 9 years",
      "10 to 14 years",
      "15 to 19 years",
      "20 to 24 years",
      "25 to 29 years",
      "30 to 34 years",
      "35 to 39 years",
      "40 to 44 years",
      "45 to 49 years",
      "50 to 54 years",
      "55 to 59 years",
      "60 to 64 years",
      "65 to 69 years",
      "70 to 74 years",
      "75 to 79 years",
      "80 to 84 years",
      "85 to 89 years",
      "90 years and over"
    ),
    weights = c(
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
    ),
    metadata = list(
      last_updated = "June 2026"
    )
  )

S7::method(as.data.frame, standard_population) <- function(x, ...) {
  rlang::check_dots_empty()
  data.frame(strata = x@strata, weights = x@weights)
}
