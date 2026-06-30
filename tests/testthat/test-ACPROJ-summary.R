test_that("acproj summary requires object class acproj", {
  object <- list(
    predictions = matrix(1, 5, 5),
    pyr = matrix(1, 5, 5)
  )

  expect_error(
    summary.acproj(object),
    "Variable \"object\" must be of type \"acproj\"",
    fixed = TRUE
  )
})

test_that("acproj summary works with printpred", {
  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    cuttrd = 0.04,
    shortp = 0,
    cuttrend = 0,
    gofpvalue = 0.9,
    distribution = "Poisson",
    startage = 1,
    glm = list("a", "b", "c")
  )
  class(acproj_obj) <- "acproj"
  attr(acproj_obj, "Call") <- as.call(str2lang("acproj(cdat, pdat)"))

  out <- capture.output(summary.acproj(acproj_obj, printpred = TRUE))

  expect_match(
    out,
    "Observed and predicted values:\\(observations up to \\)",
    all = FALSE
  )

  expect_length(grep("\\[[0-9]+,\\]\\s+36\\s+60\\s+86\\s+110", out), 19)
})

test_that("acproj summary default produces correct output", {
  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    cuttrd = 0.04,
    shortp = 0,
    cuttrend = 0,
    gofpvalue = 0.9,
    distribution = "Poisson",
    startage = 1,
    glm = list("a", "b", "c")
  )
  class(acproj_obj) <- "acproj"
  attr(acproj_obj, "Call") <- as.call(str2lang("acproj(cdat, pdat)"))

  out <- capture.output(summary.acproj(acproj_obj, printcall = TRUE))

  expect_match(out, "Prediction done with:", all = FALSE)
  expect_match(out, "Method:\\s+Age-Cohort Model", all = FALSE)
  expect_match(
    out,
    "Number of periods predicted \\(nopred\\):\\s+1",
    all = FALSE
  )
  expect_match(
    out,
    "Trend used in new cohort estimation \\(cuttrend\\):\\s+0",
    all = FALSE
  )
  expect_match(
    out,
    "Number of periods used in estimate \\(noperiod\\):\\s+3",
    all = FALSE
  )
  expect_match(
    out,
    "Distribution function of regression:\\s+Poisson",
    all = FALSE
  )
  expect_match(out, "P-value for goodness of fit:\\s+0.9", all = FALSE)
  expect_match(
    out,
    "First age group estimated \\(startage\\):\\s+1",
    all = FALSE
  )
})

test_that("acproj summary can create printcall", {
  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    cuttrd = 0.04,
    shortp = 0,
    cuttrend = 0,
    gofpvalue = 0.9,
    distribution = "Poisson",
    startage = 1,
    glm = list("a", "b", "c")
  )
  class(acproj_obj) <- "acproj"
  attr(acproj_obj, "Call") <- as.call(str2lang("acproj(cdat, pdat)"))

  out <- capture.output(summary.acproj(acproj_obj, printcall = TRUE))

  expect_match(out, "Call:\\s+acproj\\(cdat, pdat\\)", all = FALSE)
})

test_that("acproj summary can specify number of digits to print", {
  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    cuttrd = 0.04,
    shortp = 0,
    cuttrend = 0,
    gofpvalue = 0.9,
    distribution = "Poisson",
    startage = 1,
    glm = list("a", "b", "c")
  )
  class(acproj_obj) <- "acproj"
  attr(acproj_obj, "Call") <- as.call(str2lang("acproj(cdat, pdat)"))

  out <- capture.output(summary.acproj(
    acproj_obj,
    printpred = TRUE,
    digits = 1
  ))

  expect_length(grep("\\[[0-9]+,\\]\\s+35.5\\s+60.5\\s+85.5\\s+110.5", out), 19)
})
