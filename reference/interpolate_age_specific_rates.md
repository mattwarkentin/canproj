# Linear interpolation

Convert to annual age-specific rates by linear interpolation

## Usage

``` r
interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
