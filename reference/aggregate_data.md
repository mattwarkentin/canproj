# Data aggregation

Aggregation of data by years

## Usage

``` r
aggregate_data(cdat, pdat, nagg)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
