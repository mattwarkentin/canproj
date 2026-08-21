# hybdproj_get_predictions

Extract projection results by n-year period

## Usage

``` r
hybdproj_get_predictions(
  hybdproj.object,
  incidence = T,
  standpop = NULL,
  excludeobs = F,
  byage,
  agegroups = "all"
)
```

## Arguments

- hybdproj.object:

  An object based on the 'hybdproj()' function.

- incidence:

  Whether to give rates (`T`) or numbers (`F`).

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age group in a standard population.

- excludeobs:

  Whether to include observed values (`T` or `F`).

- byage:

  Report numbers by age groups (`T`), or use age-standardized rates
  (`F`).

- agegroups:

  Age groups to include. "all" (default) includes all age groups, and
  individual groups can be selected by group number. E.g. c(1:3, 7)
  includes the first 3 groups and the seventh group.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
