# HYBDPROJ

R functions for projection of cancer incidence/mortality using the
modified Hybrid methods. Modified Hybrid by adding choice of age-model,
cut-trend parameter, and power 5 link function.

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
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

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

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
