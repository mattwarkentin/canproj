# hybdproj_predict

Extrapolate estimated trend from the final model

## Usage

``` r
hybdproj_predict(hybdproj_estimate.object, cuttrd = 0.04, shortp = 0)
```

## Arguments

- hybdproj_estimate.object:

  An object based on the `hybdproj_estimate` function.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
