# Convert to age-standardized

Calculate age-standardized annual rates and counts

## Usage

``` r
standardize_annual_rates(rates, pdat, startp, stdpop, check = TRUE)
```

## Arguments

- rates:

  Observed and projected age-specific rates.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- startp:

  The start calendar year of projection (e.g., 2009).

- check:

  Whether to check APC for extreme differences
