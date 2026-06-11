test_that("conversion to age-standardized rates and cases works", {
  rr <- matrix(1:57, nrow = 19, ncol = 3)
  pdat <- matrix(1000 + 10 * (1:57), nrow = 19, ncol = 3)
  stdpop <- rep(0.2, 19)

  expected_matrix <- cbind(c(38, 110.2, 182.4), c(2, 7, 14))
  colnames(expected_matrix) <- c("asr", "case")

  expect_equal(asry(rr, pdat, standpop = stdpop), expected_matrix)
})

test_that("conversion protects against negative values", {
  rr <- matrix(c(rep(-1, 19), rep(0, 19), rep(1, 19)), nrow = 19, ncol = 3)
  pdat <- matrix(10000, nrow = 19, ncol = 3)

  expected_matrix <- cbind(c(0, 0, 1), c(0, 0, 2))
  colnames(expected_matrix) <- c("asr", "case")

  expect_equal(asry(rr, pdat, stdpop_Canada_2021()), expected_matrix)
})
