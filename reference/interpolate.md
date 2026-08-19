# Interpolation

Calculate rates by linear interpolation

## Usage

``` r
interpolate(
  nagg,
  matrix,
  rates,
  next_rate,
  last_rate,
  npred,
  nagg_proj_year,
  ny_proj
)
```

## Arguments

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

- matrix:

  Base matrix to work off of.

- rates:

  Matrix of rates for predicted periods.

- next_rate:

  Calculated rate for extra years.

- last_rate:

  Calculated rate for extra years.

- npred:

  Number of predicted periods after aggregation.

- nagg_proj_year:

  Number of aggregated projection years.

- ny_proj:

  Number of years for preojection, not aggregated.

## Value

A [`matrix()`](https://rdrr.io/r/base/matrix.html).
