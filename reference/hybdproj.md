# HYBDPROJ: Hybrid Projections

R functions for projection of cancer incidence/mortality using the
modified Hybrid methods. Modified Hybrid by adding choice of age-model,
cut-trend parameter, and power 5 link function.

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
extracts projection results from a `hybdproj()` object.

## Usage

``` r
hybdproj(
  cdat,
  pdat,
  standpop,
  projfor = "incidence",
  nagg = NULL,
  ncase = NULL,
  cuttrd = 0.04,
  shortp = 0,
  linkfunc = "power5",
  pD = 0.05,
  pGOF = 0.05
)

# S3 method for class 'hybdproj'
get_projections(
  object,
  ...,
  cdat,
  pdat,
  startp,
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

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age group in a standard population.

- projfor:

  Specify `"incidence"` or `"mortality"` if you want ASR as a criteria
  for `nagg`.

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

- ncase:

  Minimum number of cancer cases/deaths per year for splitting data.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

- pD:

  Trend selecting criteria of p-value of drift (linear trend) term.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

- object:

  Output from `hybdproj()`.

- ...:

  Not currently used.

- startp:

  The start calendar year of projection (e.g., 2009).

- ave5:

  `ave5 = TRUE` invokes the 5-year average method when age-only model is
  selected.

- sum5:

  When the 5-year average method is used, `sum5 = TRUE` call the 5-year
  period based rate, otherwise (`sum5 = FALSE`), average the 5 rates in
  the 5 years for each age group.

## Value

`hybdproj()` returns a `list`.

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
returns a `data.frame`.
