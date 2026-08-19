# Select hybrid model

Select projection model for hybrid method

## Usage

``` r
select_hybd_method(casesB, apdata, nagg, pGOF, dnoperiods, pD)
```

## Arguments

- casesB:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- apdata:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- nagg:

  Number of years for data aggregation (by years). Default is 1-annual
  data.

- pGOF:

  Model selection criteria of p-value of goodness-of-fit.

- dnoperiods:

  Number of periods.

- pD:

  Trend selecting criteria of p-value of drift (linear trend) term.
