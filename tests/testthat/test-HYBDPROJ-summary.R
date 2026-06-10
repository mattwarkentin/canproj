test_that("hybdproj summary requires hybdproj object", {
  object <- list(
    glm = list(1:3),
    cases = matrix(5, 5, 5),
    pyr = matrix(5, 5, 5)
  )

  expect_error(
    summary.hybdproj(object),
    "Variable \"object\" must be of type \"hybdproj\""
  )
})

test_that("hybdproj summary works with defaults", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    nototper = 4,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "common-trend",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  out <- capture.output(summary.hybdproj(hybd_obj))

  expect_match(out, "Prediction done with:", all = FALSE)
  expect_match(out, "Method:\\s+Hybrid approach", all = FALSE)
  expect_match(out, "Number of prediction years:\\s+4", all = FALSE)
  expect_match(out, "First period cutting trend:\\s+0.02", all = FALSE)
  expect_match(out, "Degenerating trend per year:\\s+0.05", all = FALSE)
  expect_match(out, "Projection base \\(years\\):\\s+13", all = FALSE)
  expect_match(out, "Aggregating years \\(nagg\\):\\s+4", all = FALSE)
  expect_match(out, "Age-cases per year \\(ncase\\):\\s+1", all = FALSE)
  expect_match(out, "Model for regression:\\s+common-trend", all = FALSE)
  expect_match(out, "Link function for GLM:\\s+power5", all = FALSE)
  expect_match(out, "P-value for goodness of fit:\\s+0.8123", all = FALSE)
  expect_match(out, "Age group for regression:\\s+1,2,3,4,5", all = FALSE)
  expect_match(out, "Age group for average method:\\s+1,2,3,6,7,8", all = FALSE)
})

test_that("hybdproj summary works with printcall", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    nototper = 4,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "common-trend",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"
  attr(hybd_obj, "Call") <- as.call(str2lang("hybdproj(cdat, pdat)"))

  out <- capture.output(summary.hybdproj(hybd_obj, printcall = TRUE))

  expect_match(out, "Call:\\s+hybdproj\\(cdat, pdat\\)", all = FALSE)
})

test_that("hybdproj summary can print predicted case numbers", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    nototper = 4,
    noobsper = 2,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "common-trend",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  out <- capture.output(summary.hybdproj(hybd_obj, printpred = TRUE))

  expect_match(
    out,
    "Predicted number of cases or deaths:\\(observations up to \\)",
    all = FALSE
  )
  expect_length(grep("\\[[0-9]+,\\]\\s+86\\s+110", out), 19)
})

test_that("hybdproj summary can adjust decimal places displayed", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35.5, nrow = 19, ncol = 4),
    nototper = 4,
    noobsper = 2,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "common-trend",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  out <- capture.output(summary.hybdproj(
    hybd_obj,
    printpred = TRUE,
    digits = 1
  ))

  expect_match(
    out,
    "Predicted number of cases or deaths:\\(observations up to \\)",
    all = FALSE
  )
  expect_length(grep("\\[[0-9]+,\\]\\s+85.5\\s+110.5", out), 19)
})
