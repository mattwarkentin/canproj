# hybdproj_estimate

Define projection base, model selection, and fitting

## Usage

``` r
hybdproj_estimate(
  cases,
  pyr,
  nagg,
  ncase,
  linkfunc = "power5",
  pD = 0.05,
  pGOF = 0.05
)
```

## Arguments

- cases:

  `data.frame` with number of cases in `nagg`-year period by ascending
  age groups in row.

- pyr:

  `data.frame` with observed and projected population size in
  `nagg`-year. period by ascending age groups in row.

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

- ncase:

  Minimum number of cancer cases/deaths per year for splitting data.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

- pD:

  Trend selecting criteria of p-value of drift (linear trend) term.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
