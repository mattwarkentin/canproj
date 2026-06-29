test_that("conversion to age-standardized rates and cases works", {
  rr <- matrix(1:57, nrow = 19, ncol = 3)
  pdat <- matrix(1000 + 10 * (1:57), nrow = 19, ncol = 3)
  stdpop <- StandardPopulation(
    "dummy",
    as.character(1:19),
    c(0.1, rep(0.05, 18))
  )

  expected_matrix <- cbind(c(9.55, 28.55, 47.55), c(2, 7, 14))
  colnames(expected_matrix) <- c("asr", "case")

  expect_equal(
    standardize_annual_rates(rr, pdat, stdpop = stdpop),
    expected_matrix
  )
})

test_that("conversion protects against negative values", {
  rr <- matrix(c(rep(-1, 19), rep(0, 19), rep(1, 19)), nrow = 19, ncol = 3)
  pdat <- matrix(10000, nrow = 19, ncol = 3)

  expected_matrix <- cbind(c(0, 0, 1), c(0, 0, 2))
  colnames(expected_matrix) <- c("asr", "case")

  expect_equal(
    standardize_annual_rates(rr, pdat, stdpop_Canada_2021),
    expected_matrix
  )
})
