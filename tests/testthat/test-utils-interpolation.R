test_that("asrpy interpolates annual rates from 5-year aggregation", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100000, nrow = 19, ncol = 10)
  startp <- 2000
  nagg <- 5
  rate <- matrix(c(39:57, 134:152), nrow = 19, ncol = 2)

  out <- asrpy(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  colnames(expected) <- c(1995:2004)
  rownames(expected) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expect_equal(out, expected)
})

test_that("asrpy interpolates annual rates from 4-year aggregation", {
  cdat <- matrix(1:76, nrow = 19, ncol = 4)
  pdat <- matrix(100000, nrow = 19, ncol = 8)
  startp <- 2000
  nagg <- 4
  rate <- matrix(c(29.5:47.5, 105.5:123.5), nrow = 19, ncol = 2)

  out <- asrpy(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:152, nrow = 19, ncol = 8))
  colnames(expected) <- c(1996:2003)
  rownames(expected) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expect_equal(out, expected)
})

test_that("asrpy interpolates annual rates from 3-year aggregation", {
  cdat <- matrix(1:57, nrow = 19, ncol = 3)
  pdat <- matrix(100000, nrow = 19, ncol = 6)
  startp <- 2000
  nagg <- 3
  rate <- matrix(c(20:38, 77:95), nrow = 19, ncol = 2)

  out <- asrpy(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:114, nrow = 19, ncol = 6))
  colnames(expected) <- c(1997:2002)
  rownames(expected) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expect_equal(out, expected)
})

test_that("asrpy interpolates annual rates from 2-year aggregation", {
  cdat <- matrix(1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 2
  rate <- matrix(c(10.5:28.5, 48.5:66.5), nrow = 19, ncol = 2)

  out <- asrpy(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(1:76, nrow = 19, ncol = 4))
  colnames(expected) <- c(1998:2001)
  rownames(expected) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expect_equal(out, expected)
})

test_that("interpolation without aggregation returns original rates", {
  cdat <- matrix(2 * 1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 1
  rate <- matrix(2 * 1:76, nrow = 19, ncol = 4)

  out <- asrpy(rate, cdat, pdat, startp, nagg)

  expected <- data.frame(matrix(2 * 1:76, nrow = 19, ncol = 4))
  colnames(expected) <- c(1998:2001)
  rownames(expected) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expect_equal(out, expected)
})

test_that("nagg must be integer from 1 to 5", {
  cdat <- matrix(2 * 1:38, nrow = 19, ncol = 2)
  pdat <- matrix(100000, nrow = 19, ncol = 4)
  startp <- 2000
  nagg <- 0
  rate <- matrix(2 * 1:76, nrow = 19, ncol = 4)

  expect_error(asrpy(rate, cdat, pdat, startp, nagg))

  nagg <- 6
  expect_error(asrpy(rate, cdat, pdat, startp, nagg))

  nagg <- 3.5
  expect_error(
    asrpy(rate, cdat, pdat, startp, nagg),
    "Years of aggregation \"nagg\" must be integer between 1 and 5"
  )
})
