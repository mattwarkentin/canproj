# ACPROJ

R functions for projection of cancer incidence/mortality. Revising and
combining nordpred and Osmond's to extrapolation cohort when no drift
appears from nordpred.

## Usage

``` r
acproj(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  startage = NULL,
  cuttrd = 0.04,
  shortp = 0,
  pGOF = 0.05,
  linkfunc = "power5"
)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- projfor:

  Specify `"incidence"` or `"mortality"` if you want ASR as a criteria
  for `nagg`.

- n5case:

  Minimum number of cancer cases/deaths per 5 years for splitting data.

- startage:

  Youngest age group to include in the GLM. Default (`NULL`) picks based
  on mean cases.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
