# acproj_predict

Extrapolate estimated trend from age-cohort model

## Usage

``` r
acproj_predict(acproj_estimate.object, cuttrd = 0.04, shortp = 0)
```

## Arguments

- acproj_estimate.object:

  An object based on the
  [`acproj_estimate()`](https://mattwarkentin.github.io/canproj/reference/acproj_estimate.md)
  function.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
