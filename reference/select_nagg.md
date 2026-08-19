# Select nagg

Select number of years to aggregate data by.

## Usage

``` r
select_nagg(cdat, pdat, standpop, projfor)
```

## Arguments

- cdat:

  (age groups) \* N (years) historical cancer data, 15\<=N\<=125.

- pdat:

  (age groups) \* (N + M) (years) observed and projected population,
  5\<=M\<=25.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age groups in a standard population.

- projfor:

  Specify `"incidence"` or `"mortality"` if you want ASR as a criteria
  for `nagg`.

## Value

A whole number.
