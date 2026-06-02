test_that("adpcproj estimate must have 19 age groups", {
  cases <- matrix(5, nrow = 12, ncol = 3)
  pyr <- matrix(10, nrow = 34, ncol = 4)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 3,
      startestage = 1,
      pGOF = 1,
      linkfunc = "power5"
    ),
    "\"cases\" and \"pyr\" must have data for 19 age groups",
    fixed = TRUE
  )
})

test_that("adpcproj estimate must have population data for all periods in cases", {
  cases <- matrix(5, nrow = 19, ncol = 4)
  pyr <- matrix(15, nrow = 19, ncol = 3)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 4,
      startestage = 1,
      pGOF = 1,
      linkfunc = "power5"
    ),
    "\"pyr\" must include information about all periods in \"cases\""
  )
})

test_that("adpcproj estimate must have population data for projection", {
  cases <- matrix(5, nrow = 19, ncol = 4)
  pyr <- matrix(15, nrow = 19, ncol = 4)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 4,
      startestage = 1,
      pGOF = 1,
      linkfunc = "power5"
    ),
    "\"pyr\" must include information on future rates"
  )
})

test_that("adpcproj estimate can't project more than 6 periods", {
  cases <- matrix(5, nrow = 19, ncol = 2)
  pyr <- matrix(15, nrow = 19, ncol = 9)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 4,
      startestage = 1,
      pGOF = 1,
      linkfunc = "power5"
    ),
    "Package can not project more than 6 periods"
  )
})

test_that("adpcproj estimate requires at least 3 periods", {
  cases <- matrix(5, nrow = 19, ncol = 3)
  pyr <- matrix(15, nrow = 19, ncol = 4)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 2,
      startestage = 1,
      pGOF = 1,
      linkfunc = "power5"
    ),
    "\"noperiod\" must be 3 or larger"
  )
})

test_that("adpcproj estimate works with power5 link function (default)", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- adpcproj.estimate(cases, pyr, noperiod = 3, startestage = 1, pGOF = 0)

  expect_s3_class(out, "adpcproj.estimate")
  expect_equal(out$cases, cases)
  expect_equal(out$pyr, pyr)
  expect_equal(out$noperiod, 3)
  expect_equal(out$linkfunc, "power5")
  expect_equal(out$startestage, 1)
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$gofpvalue, 0.9917539, tolerance = 0.000001)
  expect_equal(out$suggestionrecent, TRUE)
  expect_equal(out$pvaluerecent, 1.497779e-07, tolerance = 1e-08)

  expect_s3_class(out$glm, "glm")
  expect_equal(
    unname(out$glm$coefficients[16]),
    0.2658555,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with log link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "log"
  )

  expect_equal(out$linkfunc, "log")
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$gofpvalue, 0.8797211, tolerance = 0.000001)

  expect_equal(unname(out$glm$coefficients[6]), -6.420765, tolerance = 0.000001)
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with sqrt link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- suppressWarnings(adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "sqrt"
  ))

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$gofpvalue, 1)

  expect_equal(unname(out$glm$coefficients[12]), 0.036293, tolerance = 0.00001)
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with identity link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- suppressWarnings(adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 0,
    linkfunc = "identity"
  ))

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$gofpvalue, 1)

  expect_equal(
    unname(out$glm$coefficients[3]),
    -0.001488186,
    tolerance = 0.00001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("adpcproj estimate fails with unknown link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  expect_error(
    adpcproj.estimate(
      cases,
      pyr,
      noperiod = 3,
      startestage = 1,
      pGOF = 0,
      linkfunc = "function"
    ),
    "Unknown \"linkfunc\"",
    fixed = TRUE
  )
})

test_that("adpcproj estimate works with negative-binomial glm and power5 link function", {
  cases <- matrix(floor(4:60 / 3) * floor(5:23 / 5) + 10, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 1,
    linkfunc = "power5"
  )

  expect_equal(out$linkfunc, "power5")
  expect_equal(out$distribution, "Negative-binomial")
  expect_equal(out$gofpvalue, 1, tolerance = 0.001)

  expect_equal(
    unname(out$glm$coefficients[17]),
    0.2906399,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with negative-binomial glm and log link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 1,
    linkfunc = "log"
  )

  expect_equal(out$linkfunc, "log")
  expect_equal(out$distribution, "Negative-binomial")
  expect_equal(out$gofpvalue, 0.8797287, tolerance = 0.000001)

  expect_equal(unname(out$glm$coefficients[2]), -7.662173, tolerance = 0.000001)
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with negative-binomial glm and sqrt link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- suppressWarnings(adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 2,
    linkfunc = "sqrt"
  ))

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$distribution, "Negative-binomial")
  expect_equal(out$gofpvalue, 1, tolerance = 0.00001)

  expect_equal(
    unname(out$glm$coefficients[12]),
    0.03629267,
    tolerance = 0.00001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("adpcproj estimate works with negative-binomial glm and identity link function", {
  cases <- matrix(1:57, nrow = 19, ncol = 3)
  pyr <- matrix(10000 + 5 * (1:76), nrow = 19, ncol = 4)

  out <- suppressWarnings(adpcproj.estimate(
    cases,
    pyr,
    noperiod = 3,
    startestage = 1,
    pGOF = 2,
    linkfunc = "identity"
  ))

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$distribution, "Negative-binomial")
  expect_equal(out$gofpvalue, 1, tolerance = 0.00001)

  expect_equal(
    unname(out$glm$coefficients[15]),
    -0.0002991863,
    tolerance = 0.00001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})
