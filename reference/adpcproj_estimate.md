# adpcproj_estimate

Fit age-drift-period-cohort models

## Usage

``` r
adpcproj_estimate(
  cases,
  pyr,
  noperiod,
  startage,
  pGOF = 0.05,
  linkfunc = "power5"
)
```

## Arguments

- cases:

  `data.frame` with number of cases in `nagg`-year period by ascending
  age groups in row.

- pyr:

  `data.frame` with observed and projected population size in
  `nagg`-year.

- noperiod:

  Number of 5-year periods in historical data.

- startage:

  Youngest age group to include in the GLM. Default (`NULL`) picks based
  on mean cases.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
