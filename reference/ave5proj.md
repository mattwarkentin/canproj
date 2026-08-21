# Average5: Five-Year Average Projections

R functions for projection of cancer incidence/mortality using the
average methods based on the recent 5 years data: (i) average numbers
and population sizes then calculate rates (default), or (ii) average the
calculated yearly rates (sum5=T)

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
extracts the projection results from an `ave5proj()` object.

## Usage

``` r
ave5proj(cdat, pdat, startp, sum5 = TRUE)

# S3 method for class 'ave5proj'
get_projections(object, ..., standpop = NULL)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

- sum5:

  When the 5-year average method is used, `sum5 = TRUE` call the 5-year
  period based rate, otherwise (`sum5 = FALSE`), average the 5 rates in
  the 5 years for each age group.

- object:

  Output object from `ave5proj()`.

- ...:

  Not currently used.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

## Value

`ave5proj()` returns a `list`.

[`get_projections()`](https://mattwarkentin.github.io/canproj/reference/get_projections.md)
returns a `data.frame`.
