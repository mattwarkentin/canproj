# Create cuttrend for hybdproj

Create cuttrend for hybdproj

## Usage

``` r
get_hybd_cuttrend(shortp, nonewpred, nagg, cuttrd)
```

## Arguments

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- nonewpred:

  Number of predicted periods

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.
