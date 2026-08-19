# adpcproj_get_projections

Extract annual projection results

## Usage

``` r
adpcproj_get_projections(cdat, pdat, startp, adpcproj.object, standpop = NULL)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

- adpcproj.object:

  An object based on the
  [`adpcproj()`](https://mattwarkentin.github.io/canproj/reference/adpcproj.md)
  function.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
