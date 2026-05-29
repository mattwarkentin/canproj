test_that("acproj prediction requires class type acproj.estimate", {
  object <- list(
    glm = list(1:3),
    cases = matrix(5, 5, 5),
    pyr = matrix(5, 5, 5)
  )

  expect_error(
    acproj.prediction(object),
    "Variable \"acproj.estimate.object\" must be of type \"acproj.estimate\"",
    fixed = TRUE
  )
})

test_that("acproj prediction works with identity link function (default)", {
  glm <- list(
    coefficients = c(
      seq(from = 0.06, to = 0.24, by = 0.01),
      seq(from = -0.1, to = 0.09, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    glm = glm,
    cases = data.frame(matrix(1:19 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    maxc = 21,
    midc = 11,
    noperiod = 3,
    linkfunc = "identity",
    startestage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8
  )
  class(estimate_obj) <- "acproj.estimate"

  out <- acproj.prediction(estimate_obj, cuttrd = 0.1, shortp = 0.01)

  expect_equal(out$glm, estimate_obj$glm)
  expect_equal(out$pyr, estimate_obj$pyr)
  expect_equal(out$nopred, 1)
  expect_equal(out$noperiod, 3)
  expect_equal(out$cuttrd, 0.1)
  expect_equal(out$shortp, 0.01)
  expect_equal(out$gofpvalue, 0.8)
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$startestage, 1)

  expect_equal(
    out$predictions["45-49", "X4"],
    1652.8000,
    tolerance = 0.0001
  )
  expect_snapshot_value(
    out$predictions,
    style = "json2",
    tolerance = 0.0001
  )
})

test_that("acproj prediction works with log link function", {
  glm <- list(
    coefficients = c(
      seq(from = 0.06, to = 0.24, by = 0.01),
      seq(from = -0.1, to = 0.09, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    glm = glm,
    cases = data.frame(matrix(1:19 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    maxc = 21,
    midc = 11,
    noperiod = 3,
    linkfunc = "log",
    startestage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8
  )
  class(estimate_obj) <- "acproj.estimate"

  out <- acproj.prediction(estimate_obj, cuttrd = 0.1, shortp = 0.01)

  expect_equal(out$predictions["80-84", "X4"], 12285.68, tolerance = 0.01)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.01)
})

test_that("acproj prediction works with sqrt link fucntion", {
  glm <- list(
    coefficients = c(
      seq(from = 0.06, to = 0.24, by = 0.01),
      seq(from = -0.1, to = 0.09, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    glm = glm,
    cases = data.frame(matrix(1:19 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    maxc = 21,
    midc = 11,
    noperiod = 3,
    linkfunc = "sqrt",
    startestage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8
  )
  class(estimate_obj) <- "acproj.estimate"

  out <- acproj.prediction(estimate_obj, cuttrd = 0.1, shortp = 0.01)

  expect_equal(out$predictions["30-34", "X4"], 264.0640, tolerance = 0.0001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.0001)
})

test_that("acproj prediction works with power5 link function", {
  glm <- list(
    coefficients = c(
      seq(from = 0.06, to = 0.24, by = 0.01),
      seq(from = -0.1, to = 0.09, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    glm = glm,
    cases = data.frame(matrix(1:19 + 20, nrow = 19, ncol = 3)),
    pyr = data.frame(matrix(
      seq(from = 10000, to = 10375, by = 5),
      nrow = 19,
      ncol = 4
    )),
    maxc = 21,
    midc = 11,
    noperiod = 3,
    linkfunc = "power5",
    startestage = 1,
    distribution = "Poisson",
    gofpvalue = 0.8
  )
  class(estimate_obj) <- "acproj.estimate"

  out <- acproj.prediction(estimate_obj, cuttrd = 0.1, shortp = 0.01)

  expect_equal(
    out$predictions["60-64", "X4"],
    1.46884207,
    tolerance = 0.00000001
  )
  expect_snapshot_value(
    out$predictions,
    style = "json2",
    tolerance = 0.00000001
  )
})
