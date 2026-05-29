test_that("acproj estimate works with power5 link function (default)", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(cases, pyr, noperiod = 3, startestage = 1, pGOF = 0)

  expect_equal(out$cases, cases)
  expect_equal(out$pyr, pyr)
  expect_equal(out$maxc, 21)
  expect_equal(out$midc, 11)
  expect_equal(out$noperiod, 3)
  expect_equal(out$linkfunc, "power5")
  expect_equal(out$startestage, 1)
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$gofpvalue, 0)

  expect_equal(
    unname(out$glm$coefficients[5]),
    0.04142698,
    tolerance = 0.0000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj estimate works with log link function", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "log"
  )

  expect_equal(out$linkfunc, "log")
  expect_equal(out$distribution, "Poisson")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -8.896974,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj estimate works with sqrt link function", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "sqrt"
  )

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$distribution, "Poisson")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -252.5895,
    tolerance = 0.0001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.0001
  )
})

test_that("acproj estimate works with identity link function", {
  cases <- matrix(1:57 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "identity"
  )

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$distribution, "Poisson")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -0.03531735,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.0000001
  )
})

test_that("acproj estimate fails with unknown link function", {
  cases <- matrix(1:57 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  expect_error(
    acproj.estimate(
      cases,
      pyr,
      noperiod = 3,
      startestage = 1,
      pGOF = 0,
      linkfunc = "other"
    ),
    "Unknown \"linkfunc\"",
    fixed = TRUE
  )
})

test_that("acproj estimate works with negative-binomial glm and power5 link function", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(cases, pyr, noperiod = 3, startestage = 1, pGOF = 1)

  expect_equal(out$linkfunc, "power5")
  expect_equal(out$distribution, "Negative-binomial")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -0.01301036,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj estimate works with negative-binomial glm and log link function", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 1,
    linkfunc = "log"
  )

  expect_equal(out$linkfunc, "log")
  expect_equal(out$distribution, "Negative-binomial")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -9.731824,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj estimate works with negative-binomial glm and sqrt link function", {
  cases <- matrix(0:56 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 1,
    linkfunc = "sqrt"
  )

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$distribution, "Negative-binomial")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -0.1325743,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj estimate works with negative-binomial glm and identity link function", {
  cases <- matrix(1:57 * 25, nrow = 19, ncol = 3)
  pyr <- matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4)

  out <- acproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 10,
    linkfunc = "identity"
  )

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$distribution, "Negative-binomial")

  expect_equal(
    unname(out$glm$coefficients[5]),
    -0.03531739,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.0000001
  )
})
