test_that("can interpolate age-specific rates from 5-year aggregation", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100000, nrow = 19, ncol = 10)
  startp <- 2000
  nagg <- 5
  rate <- matrix(c(39:57, 134:152), nrow = 19, ncol = 2)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  colnames(expected) <- c(1995:2004)
  rownames(expected) <- 1:19

  expect_equal(out, expected)

  pdat <- matrix(100000, 19, 11)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[1, 11], 191)

  pdat <- matrix(100000, 19, 12)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[5, 12], 214)

  pdat <- matrix(100000, 19, 13)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[14, 13], 242)

  pdat <- matrix(100000, 19, 14)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[10, 14], 238)
})

test_that("can interpolate age-specific rates from 4-year aggregation", {
  cdat <- matrix(1:76, nrow = 19, ncol = 4)
  pdat <- matrix(100000, nrow = 19, ncol = 8)
  startp <- 2000
  nagg <- 4
  rate <- matrix(c(29.5:47.5, 105.5:123.5), nrow = 19, ncol = 2)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:152, nrow = 19, ncol = 8))
  colnames(expected) <- c(1996:2003)
  rownames(expected) <- 1:19

  expect_equal(out, expected)

  pdat <- matrix(100000, 19, 9)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[8, 9], 160)

  pdat <- matrix(100000, 19, 10)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[18, 10], 189)

  pdat <- matrix(100000, 19, 11)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[10, 11], 190.5)
})

test_that("asrpy interpolates annual rates from 3-year aggregation", {
  cdat <- matrix(1:57, nrow = 19, ncol = 3)
  pdat <- matrix(100000, nrow = 19, ncol = 6)
  startp <- 2000
  nagg <- 3
  rate <- matrix(c(20:38, 77:95), nrow = 19, ncol = 2)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:114, nrow = 19, ncol = 6))
  colnames(expected) <- c(1997:2002)
  rownames(expected) <- 1:19

  expect_equal(out, expected)

  pdat <- matrix(100000, 19, 7)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[4, 7], 118)

  pdat <- matrix(100000, 19, 8)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[15, 7], 129)
})

test_that("can interpolate age-specific rates from 2-year aggregation", {
  cdat <- matrix(1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 2
  rate <- matrix(c(10.5:28.5, 48.5:66.5), nrow = 19, ncol = 2)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:76, nrow = 19, ncol = 4))
  colnames(expected) <- c(1998:2001)
  rownames(expected) <- 1:19

  expect_equal(out, expected)

  pdat <- matrix(100000, 19, 5)
  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)
  expect_equal(out[10, 5], 86)
})

test_that("interpolation without aggregation returns original rates", {
  cdat <- matrix(2 * 1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 1
  rate <- matrix(2 * 1:76, nrow = 19, ncol = 4)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(2 * 1:76, nrow = 19, ncol = 4))
  colnames(expected) <- c(1998:2001)
  rownames(expected) <- 1:19

  expect_equal(out, expected)
})

test_that("nagg must be integer from 1 to 5", {
  cdat <- matrix(2 * 1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 0
  rate <- matrix(2 * 1:76, nrow = 19, ncol = 4)

  expect_error(interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg))

  nagg <- 6
  expect_error(interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg))

  nagg <- 3.5
  expect_error(
    interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg),
    "Years of aggregation \"nagg\" must be integer between 1 and 5"
  )
})

test_that("interpolation accomodates varying age groups", {
  cdat <- matrix(1:50, nrow = 10, ncol = 5)
  pdat <- matrix(100000, nrow = 10, ncol = 10)
  startp <- 2000
  nagg <- 5
  rate <- matrix(c(39:48, 124:133), nrow = 10, ncol = 2)

  out <- interpolate_age_specific_rates(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(
    c(1:50, 90:99, 107:116, 124:133, 141:150, 158:167),
    nrow = 10,
    ncol = 10
  ))
  colnames(expected) <- c(1995:2004)
  rownames(expected) <- 1:10

  expect_equal(out, expected)
})
