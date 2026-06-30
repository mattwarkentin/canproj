test_that("adpcproj prediction requires object class adpcproj.estimate", {
  object <- list(
    cases = matrix(5, 5, 5),
    pyr = matrix(5, 5, 5),
    glm = list(1:5)
  )

  expect_error(
    adpcproj.prediction(object, 1, TRUE),
    "Variable \"adpcproj.estimate.object\" must be of type \"adpcproj.estimate\"",
    fixed = TRUE
  )
})

test_that("adpcproj prediction startuseage must greater or equal to start stage", {
  adpc_est <- list(
    cases = matrix(5, 5, 5),
    pyr = matrix(5, 5, 5),
    glm = list(1:5),
    startage = 5
  )
  class(adpc_est) <- "adpcproj.estimate"

  expect_error(
    adpcproj.prediction(adpc_est, startuseage = 2, TRUE),
    "\"startuseage\" is set too low compared to \"startage\"",
    fixed = TRUE
  )
})

test_that("adpcproj prediction works with recent trend and power5 link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "power5",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = TRUE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = TRUE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$recent, TRUE)

  expect_equal(out$predictions["3", "X4"], 1.948014, tolerance = 0.00001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("adpcproj prediction works with recent trend and log link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "log",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = TRUE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = TRUE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$recent, TRUE)

  expect_equal(out$predictions["8", "X4"], 12668.69, tolerance = 0.01)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.01)
})

test_that("adpcproj prediction works with recent trend and sqrt link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "sqrt",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = TRUE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = TRUE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$recent, TRUE)

  expect_equal(out$predictions["16", "X4"], 622.1128, tolerance = 0.0001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.0001)
})

test_that("adpcproj prediction works with recent trend and identity link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "identity",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = TRUE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = TRUE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$recent, TRUE)

  expect_equal(out$predictions["1", "X4"], 1748.964, tolerance = 0.001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.001)
})

test_that("adpcproj prediction works with power5 link function", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "power5",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = FALSE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = FALSE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$predictions["9", "X4"], 5.940696, tolerance = 0.0000001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.0000001)
})

test_that("adpcproj prediction works with log link function", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "log",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = FALSE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = FALSE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$predictions["18", "X4"], 13582.97, tolerance = 0.01)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.01)
})

test_that("adpcproj prediction works with sqrt link function", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "sqrt",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = FALSE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = FALSE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(out$predictions["14", "X4"], 646.3576, tolerance = 0.0001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.0001)
})

test_that("adpcproj prediction works with identity link function (default)", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.125, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:19),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:21)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:57 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "identity",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = FALSE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = FALSE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_s3_class(out, "adpcproj")
  expect_equal(out$glm, estimate_adpc$glm)
  expect_equal(out$pyr, estimate_adpc$pyr)

  expect_equal(out$nopred, 1)
  expect_equal(out$noperiod, 3)
  expect_equal(out$gofpvalue, 0.8)
  expect_equal(out$recent, FALSE)
  expect_equal(out$pvaluerecent, 0.01)
  expect_equal(out$cuttrd, 0.05)
  expect_equal(out$shortp, 0.01)
  expect_equal(out$cuttrend, 0.01)
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$startuseage, 1)
  expect_equal(out$startage, 1)

  expect_equal(out$predictions["13", "X4"], 2533.490, tolerance = 0.001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.001)
})

test_that("adpcproj prediction works with varying age groups", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.11, by = 0.01),
      seq(from = 0.01, to = 0.08, by = 0.005)
    ),
    method = "fake"
  )
  attributes(glm$coefficients)$names <- c(
    paste0("as.factor(Age)", 1:10),
    "Period",
    paste0("as.factor(Period)", 2:3),
    paste0("as.factor(Cohort)", 2:12)
  )

  estimate_adpc <- list(
    glm = glm,
    cases = data.frame(matrix(1:30 + 20, nrow = 10, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10195, by = 5),
      nrow = 10,
      ncol = 4
    )),
    noperiod = 3,
    linkfunc = "identity",
    startage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8,
    suggestionrecent = FALSE,
    pvaluerecent = 0.01
  )
  class(estimate_adpc) <- "adpcproj.estimate"

  out <- adpcproj.prediction(
    estimate_adpc,
    startuseage = 1,
    recent = FALSE,
    shortp = 0.01,
    cuttrd = 0.05
  )

  expect_equal(nrow(out$predictions), 10)
  expect_equal(out$predictions["7", "X4"], 1729.582, tolerance = 0.001)
})
