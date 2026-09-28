# Project Using Constant Rates

Project future cancer cases by applying age-specific cancer rates from
the most recent year of observed data to projected future population
distributions. The function assumes that the most recent observed cancer
rates remain constant over the projection period and applies these rates
to the corresponding age strata in each future year.

## Usage

``` r
project_constant(cdat, pdat, startp)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
