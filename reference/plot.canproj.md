# plot.canproj

Create the graph of the observed and projected age-standardized rates
from a canproj object.

## Usage

``` r
# S3 method for class 'canproj'
plot(
  x,
  standpop,
  startplot = 1,
  xlab = "Calendar Year",
  ylab = "Rate per 100,000 people",
  main = "",
  lty = c(1, 2),
  col = c("black", "azure4"),
  ...
)
```

## Arguments

- x:

  An object based on the 'canproj()' function.

- standpop:

  A `StandardPopulation` object that provides the weights (proportions)
  for each age group in a standard population.

- startplot:

  Start for usage of data. Default (`1`) uses all years, increasing this
  number removes the oldest years from the graph.

- xlab:

  x-axis label

- ylab:

  y-axis label

- main:

  Title for graph

- lty:

  Line type. Applies to observed rates and predicted rates,
  respectively.

- col:

  Line colour. Applies to observed rates and predicted rates,
  respectively.

- ...:

  Other parameters
