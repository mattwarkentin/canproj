test_that("adpcproj summary requires object class adpcproj", {
  object <- list(
    predictions = matrix(1, 5, 5),
    pyr = matrix(1, 5, 5)
  )

  expect_error(
    summary.adpcproj(object),
    "Variable \"object\" must be of type \"adpcproj\"",
    fixed = TRUE
  )
})

test_that("adpcproj summary works with default settings", {
  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    recent = FALSE,
    pvaluerecent = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    cuttrend = 0.02,
    gofpvalue = 1,
    distribution = "Poisson",
    startestage = 1,
    startuseage = 1,
    glm = list("a", "b", "c")
  )
  class(adpcproj_obj) <- "adpcproj"
  attr(adpcproj_obj, "Call") <- as.call(str2lang("adpcproj(cdat, pdat)"))

  out <- capture.output(summary.adpcproj(adpcproj_obj))

  expect_match(out, "Prediction done with:", all = FALSE)
  expect_match(out, "Method:\\s+Age-drift-Period-Cohort Model", all = FALSE)
  expect_match(
    out,
    "Number of periods predicted \\(nopred\\):\\s+1",
    all = FALSE
  )
  expect_match(
    out,
    "Trend used in predictions \\(cuttrend\\):\\s+0.02",
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
  expect_match(out, "P-value for goodness of fit:\\s+1", all = FALSE)
  expect_match(out, "Used recent \\(recent\\):\\s+FALSE", all = FALSE)
  expect_match(out, "P-value for recent:\\s+1", all = FALSE)
  expect_match(out, "First age group used \\(startuseage\\):\\s+1", all = FALSE)
  expect_match(
    out,
    "First age group estimated \\(startestage\\):\\s+1",
    all = FALSE
  )
})

test_that("adpcproj summary works with printpred", {
  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    recent = FALSE,
    pvaluerecent = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    cuttrend = 0.02,
    gofpvalue = 1,
    distribution = "Poisson",
    startestage = 1,
    startuseage = 1,
    glm = list("a", "b", "c")
  )
  class(adpcproj_obj) <- "adpcproj"
  attr(adpcproj_obj, "Call") <- as.call(str2lang("adpcproj(cdat, pdat)"))

  out <- capture.output(summary.adpcproj(adpcproj_obj, printpred = TRUE))

  expect_match(
    out,
    "Observed and predicted values:\\(observations up to \\)",
    all = FALSE
  )
  expect_length(grep("\\[[0-9]+,\\]\\s+36\\s+60\\s+86\\s+110", out), 19)
})

test_that("adpcproj summary works with printcall", {
  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    recent = FALSE,
    pvaluerecent = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    cuttrend = 0.02,
    gofpvalue = 1,
    distribution = "Poisson",
    startestage = 1,
    startuseage = 1,
    glm = list("a", "b", "c")
  )
  class(adpcproj_obj) <- "adpcproj"
  attr(adpcproj_obj, "Call") <- as.call(str2lang("adpcproj(cdat, pdat)"))

  out <- capture.output(summary.adpcproj(adpcproj_obj, printcall = TRUE))

  expect_match(out, "Call:\\s+adpcproj\\(cdat, pdat\\)", all = FALSE)
})

test_that("adpcproj summary can specify number of digits to print", {
  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1,
    noperiod = 3,
    recent = FALSE,
    pvaluerecent = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    cuttrend = 0.02,
    gofpvalue = 1,
    distribution = "Poisson",
    startestage = 1,
    startuseage = 1,
    glm = list("a", "b", "c")
  )
  class(adpcproj_obj) <- "adpcproj"
  attr(adpcproj_obj, "Call") <- as.call(str2lang("adpcproj(cdat, pdat)"))

  out <- capture.output(summary.adpcproj(
    adpcproj_obj,
    printpred = TRUE,
    digits = 1
  ))

  expect_match(
    out,
    "Observed and predicted values:\\(observations up to \\)",
    all = FALSE
  )
  expect_length(grep("\\[[0-9]+,\\]\\s+35.5\\s+60.5\\s+85.5\\s+110.5", out), 19)
})
