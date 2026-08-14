test_that("conversion to age-standardized rates and cases works", {
  set.seed(89135)
  rr <- matrix(rnorm(190, mean = 20, 1), nrow = 19, ncol = 10)
  pdat <- matrix(1000 + (1:190), nrow = 19, ncol = 10)
  colnames(pdat) <- 2020:2029

  stdpop <- StandardPopulation(
    "dummy",
    as.character(1:19),
    c(0.1, rep(0.05, 18))
  )

  expected_matrix <- cbind(
    c(
      19.912387,
      20.099375,
      20.405806,
      20.069178,
      20.239549,
      20.317642,
      20.140695,
      20.10037,
      19.794404,
      19.958169
    ),
    rep(4, 10)
  )
  colnames(expected_matrix) <- c("asr", "case")
  rownames(expected_matrix) <- 2020:2029

  expect_equal(
    standardize_annual_rates(rr, pdat, 2021, stdpop = stdpop),
    expected_matrix
  )
})

test_that("conversion protects against negative values", {
  rr <- matrix(
    rep(
      c(rep(-0.5, 19), rep(0.5, 19), rep(0.1, 19), rep(0.2, 19), rep(0.2, 19)),
      2
    ),
    nrow = 19,
    ncol = 10
  )
  pdat <- matrix(10000, nrow = 19, ncol = 10)
  colnames(pdat) <- 2020:2029

  expected_matrix <- cbind(
    c(0, 0.5, 0.1, 0.2, 0.2, 0, 0.5, 0.1, 0.2, 0.2),
    c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0)
  )
  colnames(expected_matrix) <- c("asr", "case")
  rownames(expected_matrix) <- 2020:2029

  expect_equal(
    standardize_annual_rates(rr, pdat, 2025, stdpop_Canada_2021),
    expected_matrix
  )
})
