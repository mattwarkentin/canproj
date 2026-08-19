# Standard Population

An `S7` class for creating standard populations to be used in other
functions in this package. For convenience, we provide three commonly
used standard populations for Canada (2011 and 2021) and the world (WHO
2000-2025). If a user requires a custom standard population, you must
use `StandardPopulation()` to construct the standard population and pass
the object as an argument to other functions, as needed.

## Usage

``` r
StandardPopulation(
  name = rlang::abort("@name is required."),
  strata = rlang::abort("@strata is required."),
  weights = rlang::abort("@weights is required."),
  metadata = list()
)

get_standard_population(name)

stdpop_Canada_2011

stdpop_Canada_2021

stdpop_WHO_2000_2025
```

## Arguments

- name:

  Name of the standard population.

- strata:

  Character vector of strata labels.

- weights:

  Numeric vector of strata weights.

- metadata:

  Optional. We provide this `metadata` property so that users can
  include arbitrary metadata in an object.

## Value

An `S7` `StandardPopulation` object.
