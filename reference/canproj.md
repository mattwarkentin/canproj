# Canproj: Cancer Projections

Cancer incidence and mortality projections using the Canproj method. In
brief, Canproj uses a decision-tree approach to determine the optimal
model to project cancer data into the future based on historical trends.
More information about this approach can be found here:
<https://doi.org/10.24095/hpcdp.40.9.02>.

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
extracts projection results from a `canproj()` object.

## Usage

``` r
canproj(
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
  methods = NULL,
  linkfunc = "power5",
  cuttrd = 0.04,
  shortp = 0,
  pD = 0.05,
  pGOF = 0.05
)

# S3 method for class 'canproj'
get_projections(object, ..., standpop = NULL)

# S3 method for class 'canproj'
print(x, ...)
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
  for each age group in a standard population.

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

- methods:

  User required projection method can be specified by ADPC models:
  `"nordpred"` or `"adpc-nb"`; age-cohort models: `"ac"`, `"ac-nb"`;
  age-period models: `"age-trd"`, `"com-trd"`; and average: `"ave5"`.

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

- object:

  Output object from `canproj()`.

- ...:

  Not currently used.

- x:

  A object of class `"canproj"` to print.

## Value

`canproj()` returns a named-`list` with class `"canproj"`. The `list`
contains the following:

- `annproj`: A `matrix` of age-standardized rates and case counts.

- `agsproj`: A `data.frame` of case counts for each age group and year.

- `method`: The chosen method for projection. One of `"ave5"`,
  `"nordpred"`, `"adpc-nb"`, `"ac-poi"`, `"ac-nb"`, `"age-trd-nb"`,
  `"age-trd-poi"`, `"age-only"`, or `"com-trd"`.

- `out`: The output from function call for the chosen `method` (e.g.,
  [`acproj()`](https://mattwarkentin.github.io/canproj/reference/acproj.md),
  [`adpcproj()`](https://mattwarkentin.github.io/canproj/reference/adpcproj.md),
  [`ave5proj()`](https://mattwarkentin.github.io/canproj/reference/ave5proj.md),
  [`hybdproj()`](https://mattwarkentin.github.io/canproj/reference/hybdproj.md)).

- `obsy`: Number of observed years of cancer data from `cdat`.

- `projy`: Number of projected years based on `pdat`.

- `pdPC`: A vector of p-values for the drift, period, and cohort effects
  in the `"adpc"` model.

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
returns a `data.frame`.
