# Make prediction table

Check to see if inputs to 'get prediction' are valid

## Usage

``` r
make_pred_table(object, agegroups, standpop, byage, incidence, excludeobs)
```

## Arguments

- object:

  "adpcproj", "hybdproj", or "acproj" object

- agegroups:

  Age groups to include. "all" (default) includes all age groups, and
  individual groups can be selected by group number. E.g. c(1:3, 7)
  includes the first 3 groups and the seventh group.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age group in a standard population.

- byage:

  Report numbers by age groups (`T`), or use age-standardized rates
  (`F`).

- incidence:

  Whether to give rates (`T`) or numbers (`F`).

- excludeobs:

  Whether to include observed values (`T` or `F`).
