# hybdproj_get_projections

Extract projection results.

## Usage

``` r
hybdproj_get_projections(
  cdat,
  pdat,
  startp,
  hybdproj.object,
  standpop = NULL,
  ave5 = FALSE,
  sum5 = TRUE
)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

- hybdproj.object:

  An object based on the 'hybdproj()' function.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

- ave5:

  `ave5 = TRUE` invokes the 5-year average method when age-only model is
  selected.

- sum5:

  When the 5-year average method is used, `sum5 = TRUE` call the 5-year
  period based rate, otherwise (`sum5 = FALSE`), average the 5 rates in
  the 5 years for each age group.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
