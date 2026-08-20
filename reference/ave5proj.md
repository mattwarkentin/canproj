# Average5: Five-Year Average Projections

R functions for projection of cancer incidence/mortality using the
average methods based on the recent 5 years data: (i) average numbers
and population sizes then calculate rates (default), or (ii) average the
calculated yearly rates (sum5=T)

## Usage

``` r
ave5proj(cdat, pdat, startp, sum5 = TRUE)
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

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
