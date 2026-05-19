test_that("obasr calculates annual age-standardized rates by cases", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100, nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- obasr(cdat, pdat, stdpop)

  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
  expect_equal(out[, 1], c(12800, 31800, 50800, 69800, 88800))
  expect_equal(out[, 2], c(190, 551, 912, 1273, 1634))
})

test_that("obasr calculates annual age-standardized rates by population", {
  cdat <- matrix(5, nrow = 19, ncol = 5)
  pdat <- matrix(100 + 2 * (1:95), nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- obasr(cdat, pdat, stdpop)

  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
  expect_equal(
    out[, 1],
    c(4007.581, 3068.095, 2486.432, 2090.523, 1803.528),
    tolerance = 0.001
  )
  expect_equal(out[, 2], rep(95, 5))
})

test_that("obasr uses Canada 2011 standard population for default", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100, nrow = 19, ncol = 20)

  out <- obasr(cdat, pdat)

  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
  expect_equal(
    out[, 1],
    c(8477.344, 27477.344, 46477.344, 65477.344, 84477.344),
    tolerance = 0.001
  )
  expect_equal(out[, 2], c(190, 551, 912, 1273, 1634))
})

test_that("asrsd calculates age-standardized rates and standard error", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(200 + 2 * (1:95), nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- asrsd(cdat, pdat, standpop = stdpop)

  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
  expect_equal(
    out[, 1],
    c(5584.634, 12008.157, 16806.539, 20527.776, 23498.211),
    tolerance = 0.001
  )
  expect_equal(
    out[, 2],
    c(432.9002, 565.0387, 619.3358, 642.7181, 650.9717),
    tolerance = 0.0001
  )
})

test_that("asrsd uses Canada 2011 standard population by default", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(200 + 2 * (1:95), nrow = 19, ncol = 20)

  out <- asrsd(cdat, pdat)

  expect_equal(nrow(out), 5)
  expect_equal(ncol(out), 2)
  expect_equal(
    out[, 1],
    c(3828.828, 10728.919, 15833.109, 19762.243, 22880.420),
    tolerance = 0.001
  )
  expect_equal(
    out[, 2],
    c(321.1489, 506.4917, 575.8266, 605.9329, 617.8074),
    tolerance = 0.0001
  )
})

test_that("chper calculates percentage change", {
  aspr <- matrix(0.05 + 0.001 * (1:95), nrow = 19, ncol = 5)
  pdat <- matrix(100 + 2 * (1:95), nrow = 19, ncol = 5)
  starty <- 2000
  byear <- 2001
  cyear <- 2004

  out <- chper(aspr, pdat, byear, cyear, starty = starty)

  expect_equal(
    colnames(out),
    c("ref.case", "comp.case", "overall", "risk", "p.growth", "p.aging")
  )
  expect_equal(nrow(out), 1)
  expect_equal(
    unname(unlist(out[1, ])),
    c(
      0.00238298,
      0.00703988,
      195.42337745,
      123.61664806,
      72.15189873,
      -0.34516934
    ),
    tolerance = 0.0000001
  )
})

test_that("chper works with NULL start year", {
  aspr <- data.frame(matrix(0.05 + 0.001 * (1:95), nrow = 19, ncol = 5))
  colnames(aspr) <- c(2000, 2001, 2002, 2003, 2004)

  pdat <- data.frame(matrix(100 + 2 * (1:95), nrow = 19, ncol = 5))
  colnames(pdat) <- c(2000, 2001, 2002, 2003, 2004)
  byear <- 2001
  cyear <- 2004

  out <- chper(aspr, pdat, byear, cyear)

  expect_equal(
    colnames(out),
    c("ref.case", "comp.case", "overall", "risk", "p.growth", "p.aging")
  )
  expect_equal(nrow(out), 1)
  expect_equal(
    unname(unlist(out[1, ])),
    c(
      0.00238298,
      0.00703988,
      195.42337745,
      123.61664806,
      72.15189873,
      -0.34516934
    ),
    tolerance = 0.0000001
  )
})
