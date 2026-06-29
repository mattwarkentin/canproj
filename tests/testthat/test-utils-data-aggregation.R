test_that("aggregation by 1 returns all years without aggregation", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 20)

  cases <- data.frame(
    matrix(1:19, nrow = 19, ncol = 15),
    row.names = c(1:19)
  )
  colnames(cases) <- 1:15

  pyr <- data.frame(
    matrix(1:19, nrow = 19, ncol = 20),
    row.names = c(1:19)
  )
  colnames(pyr) <- 1:20

  agg <- list(cases = cases, pyr = pyr)

  expect_equal(aggregate_data(cdat, pdat, 1), agg)
})

test_that("aggregation by a factor of nyears aggregates all years", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 20)

  cases <- data.frame(
    matrix(5 * (1:19), nrow = 19, ncol = 3),
    row.names = c(1:19)
  )
  colnames(cases) <- 1:3

  pyr <- data.frame(
    matrix(5 * (1:19), nrow = 19, ncol = 4),
    row.names = c(1:19)
  )
  colnames(pyr) <- 1:4

  agg <- list(cases = cases, pyr = pyr)

  expect_equal(aggregate_data(cdat, pdat, 5), agg)
})

test_that("aggregation by a nonfactor removes extra years at beginning", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(1:133, nrow = 19, ncol = 7)

  cases <- data.frame(
    a = 57 + 2 * (1:19),
    b = 133 + 2 * (1:19),
    row.names = c(1:19)
  )
  colnames(cases) <- 1:2

  pyr <- data.frame(
    a = 57 + 2 * (1:19),
    b = 133 + 2 * (1:19),
    c = 209 + 2 * (1:19),
    row.names = c(1:19)
  )
  colnames(pyr) <- 1:3

  agg <- list(cases = cases, pyr = pyr)

  expect_equal(aggregate_data(cdat, pdat, 2), agg)
})

test_that("population input must include projection years", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 15)

  expect_error(datagg(cdat, pdat, 1))
})

test_that("data aggregation supports other age groupings", {
  cdat <- matrix(1:10, nrow = 10, ncol = 15)
  pdat <- matrix(1000:1009, nrow = 10, ncol = 20)

  cases <- data.frame(matrix(5 * 1:10, 10, 3), row.names = 1:10)
  colnames(cases) <- 1:3

  pyr <- data.frame(matrix(5 * 0:9 + 5000, 10, 4), row.names = 1:10)
  colnames(pyr) <- 1:4

  agg <- list(cases = cases, pyr = pyr)

  expect_equal(aggregate_data(cdat, pdat, 5), agg)
})
