# Age-standardized rates and standard error

Calculate age-standardized rates and standard error

## Usage

``` r
standardize_rates(cdat, pdat, stdpop = stdpop_Canada_2021)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

## Value

A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
