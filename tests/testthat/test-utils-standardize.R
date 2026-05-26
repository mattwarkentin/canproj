test_that("obasr calculates annual age-standardized rates by cases", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100, nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- matrix(
    c(
      c(12800, 31800, 50800, 69800, 88800),
      c(190, 551, 912, 1273, 1634)
    ),
    nrow = 5,
    ncol = 2
  )
  colnames(out) <- c("asr", "case")

  expect_equal(obasr(cdat, pdat, stdpop), out)
})

test_that("obasr calculates annual age-standardized rates by population", {
  cdat <- matrix(5, nrow = 19, ncol = 5)
  pdat <- matrix(100 + 2 * (1:95), nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- matrix(
    c(
      c(4007.581, 3068.095, 2486.432, 2090.523, 1803.528),
      rep(95, 5)
    ),
    nrow = 5,
    ncol = 2
  )
  colnames(out) <- c("asr", "case")

  expect_equal(obasr(cdat, pdat, stdpop), out, tolerance = 0.001)
})

test_that("obasr uses Canada 2011 standard population for default", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(100, nrow = 19, ncol = 20)

  out <- matrix(
    c(
      c(8477.344, 27477.344, 46477.344, 65477.344, 84477.344),
      c(190, 551, 912, 1273, 1634)
    ),
    nrow = 5,
    ncol = 2
  )
  colnames(out) <- c("asr", "case")

  expect_equal(obasr(cdat, pdat), out, tolerance = 0.001)
})

test_that("asrsd calculates age-standardized rates and standard error", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(200 + 2 * (1:95), nrow = 19, ncol = 20)
  stdpop <- c(rep(0.02, 5), rep(0.04, 5), rep(0.06, 5), rep(0.1, 4))

  out <- matrix(
    c(
      c(5584.634, 12008.157, 16806.539, 20527.776, 23498.211),
      c(432.9002, 565.0387, 619.3358, 642.7181, 650.9717)
    ),
    nrow = 5,
    ncol = 2
  )
  colnames(out) <- c("asr", "asd")

  expect_equal(asrsd(cdat, pdat, standpop = stdpop), out, tolerance = 0.001)
})

test_that("asrsd uses Canada 2011 standard population by default", {
  cdat <- matrix(1:95, nrow = 19, ncol = 5)
  pdat <- matrix(200 + 2 * (1:95), nrow = 19, ncol = 20)

  out <- matrix(
    c(
      c(3828.828, 10728.919, 15833.109, 19762.243, 22880.420),
      c(321.1489, 506.4917, 575.8266, 605.9329, 617.8074)
    ),
    nrow = 5,
    ncol = 2
  )
  colnames(out) <- c("asr", "asd")

  expect_equal(asrsd(cdat, pdat), out, tolerance = 0.001)
})

test_that("chper calculates percentage change", {
  aspr <- matrix(0.05 + 0.001 * (1:95), nrow = 19, ncol = 5)
  pdat <- matrix(100 + 2 * (1:95), nrow = 19, ncol = 5)
  starty <- 2000
  byear <- 2001
  cyear <- 2004

  out <- data.frame(
    ref.case = 0.00238298,
    comp.case = 0.00703988,
    overall = 195.42337745,
    risk = 123.61664806,
    p.growth = 72.15189873,
    p.aging = -0.34516934
  )

  expect_equal(
    chper(aspr, pdat, byear, cyear, starty = starty),
    out,
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

  out <- data.frame(
    ref.case = 0.00238298,
    comp.case = 0.00703988,
    overall = 195.42337745,
    risk = 123.61664806,
    p.growth = 72.15189873,
    p.aging = -0.34516934
  )

  expect_equal(
    chper(aspr, pdat, byear, cyear),
    out,
    tolerance = 0.0000001
  )
})
