# ADPCPROJ

R functions for projection of cancer incidence/mortality. Revising
nordpred and introducing negative binomial distribution when lack of fit
appears from nordpred, additional link functions of sqrt and identity,
and settings of startage and startuseage.

## Usage

``` r
adpcproj(
  cdat,
  pdat,
  projfor = "incidence",
  n5case = NULL,
  noperiods = NULL,
  recent = NULL,
  startage = NULL,
  newcohort = FALSE,
  pGOF = 0.05,
  cuttrd = 0.04,
  shortp = 0,
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

- noperiods:

  List of candidate periods for projection base. Default (`NULL`) uses a
  goodness-of-fit test to determine if ancient periods are removed.

- recent:

  Estimate drift term from recent trend (`T`) or whole trend (`F`).
  Default (`NULL`) uses compares models to pick.

- startage:

  Youngest age group to include in the GLM. Default (`NULL`) picks based
  on mean cases.

- newcohort:

  Assign new cohort effect as `0` (`FALSE`) or the last estimated cohort
  effect (`TRUE`), default is `0`, use `TRUE` only if having evidence on
  negative new cohort effect.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
