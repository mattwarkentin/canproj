# ave5proj_get_projections

Extract the projection results

## Usage

``` r
ave5proj_get_projections(pdat, ave5proj.object, standpop = NULL)
```

## Arguments

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- ave5proj.object:

  An object based on the
  [`ave5proj()`](https://mattwarkentin.github.io/canproj/reference/ave5proj.md)
  function.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
