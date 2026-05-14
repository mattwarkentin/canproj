test_that("aggregation by 1 returns all years without aggregation", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 20)

  agg <- datagg(cdat, pdat, 1)
  expect_equal(length(agg), 2)

  expected_cases <- data.frame(
    matrix(1:19, nrow = 19, ncol = 15),
    row.names = c(1:19)
  )
  colnames(expected_cases) <- 1:15
  expect_equal(agg$cases, expected_cases)

  expected_pyr <- data.frame(
    matrix(1:19, nrow = 19, ncol = 20),
    row.names = c(1:19)
  )
  colnames(expected_pyr) <- 1:20
  expect_equal(agg$pyr, expected_pyr)
})

test_that("aggregation by a factor of nyears aggregates all years", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 20)

  agg <- datagg(cdat, pdat, 5)
  expect_equal(length(agg), 2)

  expected_cases <- data.frame(
    matrix(5 * (1:19), nrow = 19, ncol = 3),
    row.names = c(1:19)
  )
  colnames(expected_cases) <- 1:3
  expect_equal(agg$cases, expected_cases)

  expected_pyr <- data.frame(
    matrix(5 * (1:19), nrow = 19, ncol = 4),
    row.names = c(1:19)
  )
  colnames(expected_pyr) <- 1:4
  expect_equal(agg$pyr, expected_pyr)
})

test_that("aggregation by a nonfactor removes extra years at beginning", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(1:133, nrow = 19, ncol = 7)

  agg <- datagg(cdat, pdat, 2)
  expect_equal(length(agg), 2)

  expected_cases <- data.frame(
    a = 57 + 2 * (1:19),
    b = 133 + 2 * (1:19),
    row.names = c(1:19)
  )
  colnames(expected_cases) <- 1:2
  expect_equal(agg$cases, expected_cases)

  expected_pyr <- data.frame(
    a = 57 + 2 * (1:19),
    b = 133 + 2 * (1:19),
    c = 209 + 2 * (1:19),
    row.names = c(1:19)
  )
  colnames(expected_pyr) <- 1:3
  expect_equal(agg$pyr, expected_pyr)
})

test_that("population input must include projection years", {
  cdat <- matrix(1:19, nrow = 19, ncol = 15)
  pdat <- matrix(1:19, nrow = 19, ncol = 15)

  expect_error(datagg(cdat, pdat, 1))
})
