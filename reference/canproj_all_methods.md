# Fit all Canproj models

Fits all of the available
[`canproj()`](https://mattwarkentin.github.io/canproj/reference/canproj.md)
methods. This function accepts all the same arguments as
[`canproj()`](https://mattwarkentin.github.io/canproj/reference/canproj.md)
and they are passed on to each function call. The methods fit by this
model include: `"nordpred"`, `"adpc-nb"`, `"ac-poi"`, `"ac-nb"`,
`"age-trd-poi"`, `"age-trd-nb"`, `"com-trd"`, `"age-only"`, and
`"ave5"`.

## Usage

``` r
canproj_all_methods(
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

- startage:

  Youngest age group to include in the GLM. Default (`NULL`) picks based
  on mean cases.

- newcohort:

  Assign new cohort effect as `0` (`FALSE`) or the last estimated cohort
  effect (`TRUE`), default is `0`, use `TRUE` only if having evidence on
  negative new cohort effect.

- ave5:

  `ave5 = TRUE` invokes the 5-year average method when age-only model is
  selected.

- sum5:

  When the 5-year average method is used, `sum5 = TRUE` call the 5-year
  period based rate, otherwise (`sum5 = FALSE`), average the 5 rates in
  the 5 years for each age group.

- linkfunc:

  Link function. Default is `"power5"`. Can be one of `"log"`, `"sqrt"`,
  or `"identity"`.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- pD:

  Trend selecting criteria of p-value of drift (linear trend) term.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

## Value

A `list`. The `list` contains the following items:

- `<selected_method>`: Projection method selected by
  [`canproj()`](https://mattwarkentin.github.io/canproj/reference/canproj.md).

- `nordpred`: `canproj` object produced by the `"nordpred"` method.

- `adpc-nb`: `canproj` object produced by the `"adpc-nb"` method.

- `ac-poi`: `canproj` object produced by the `"ac-poi"` method.

- `ac-nb`:`canproj` object produced by the `"ac-nb"` method.

- `a-s-nb`: `canproj` object produced by the `"a-s-nb"` method.

- `a-s-poi`: `canproj` object produced by the `"a-s-poi"` method.

- `c-t`: `canproj` object produced by the `"c-t"` method.

- `average`: `canproj` object produced by the `"average"` method.

- `ave5`: `canproj` object produced by the `"ave5"` method.
