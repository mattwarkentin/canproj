# adpcproj_predict

Extrapolate estimated trend from age-drift-period-cohort model

## Usage

``` r
adpcproj_predict(
  adpcproj_estimate.object,
  startuseage,
  recent,
  shortp = 0,
  cuttrd = 0.04
)
```

## Arguments

- adpcproj_estimate.object:

  An object based on the
  [`adpcproj_estimate()`](https://mattwarkentin.github.io/canproj/reference/adpcproj_estimate.md)
  function.

- startuseage:

  Youngest age group to use estimates from the GLM for projection.

- recent:

  Indicate estimated drift term from recent trend (`T`) or whole trend
  (`F`).

- shortp:

  Attenuation percent of drift term or slope for the first 5 projection
  years.

- cuttrd:

  Degenerating percent of trends per year after 5 years (`shortp = 0`)
  or the first projection year.

## Value

A [`list()`](https://rdrr.io/r/base/list.html).
