test_that("hybdproj estimate works", {
  cases <- matrix(2 * 1:19 + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase)

  expect_equal(out$cases, cases)
  expect_equal(out$pyr, pyr)

  expect_equal(out$lastper, 5)
  expect_equal(out$cuty, 1)
  expect_equal(out$noperiod, 5)
  expect_equal(out$noyearagg, 4)
  expect_equal(out$nocaseagp, 5)
  expect_equal(out$agrpmod, 1:19)
  expect_equal(out$projbase, 20)
})

test_that("hybdproj estimate works with varying age groups", {
  cases <- matrix(2 * 1:10 + 20, nrow = 10, ncol = 5)
  pyr <- matrix(50001:50060, nrow = 10, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase)

  expect_equal(out$agrpmod, 1:10)
  expect_equal(length(out$glm$coefficients), 10)
})

test_that("hybdproj estimate works with average model and power5 link func", {
  cases <- matrix(2 * 1:19 + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "power5")

  expect_equal(out$finalmod, "average")
  expect_equal(out$linkfunc, "power5")
  expect_equal(unname(out$glm$coefficients[12]), 0.244799, tolerance = 0.000001)

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with average model and log link func", {
  cases <- matrix(2 * 1:19 + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "log")

  expect_equal(out$finalmod, "average")
  expect_equal(out$linkfunc, "log")
  expect_equal(
    unname(out$glm$coefficients[7]),
    -7.294317,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with average model and sqrt link func", {
  cases <- matrix(2 * 1:19 + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "sqrt")

  expect_equal(out$finalmod, "average")
  expect_equal(out$linkfunc, "sqrt")
  expect_equal(
    unname(out$glm$coefficients[4]),
    0.02365439,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with average model and identity link func", {
  cases <- matrix(2 * 1:19 + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "identity")

  expect_equal(out$finalmod, "average")
  expect_equal(out$linkfunc, "identity")
  expect_equal(
    unname(out$glm$coefficients[17]),
    0.001078814,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with common-trend model and power5 link func", {
  cases <- matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "power5")

  expect_equal(out$finalmod, "common-trend")
  expect_equal(out$linkfunc, "power5")
  expect_equal(
    unname(out$glm$coefficients[10]),
    0.2072827,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with common-trend model and log link func", {
  cases <- matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "log")

  expect_equal(out$finalmod, "common-trend")
  expect_equal(out$linkfunc, "log")
  expect_equal(unname(out$glm$coefficients[2]), -7.866615, tolerance = 0.000001)

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with common-trend model and sqrt link func", {
  cases <- matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "sqrt")

  expect_equal(out$finalmod, "common-trend")
  expect_equal(out$linkfunc, "sqrt")
  expect_equal(
    unname(out$glm$coefficients[16]),
    0.01953667,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with common-trend model and identity link func", {
  cases <- matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "identity")

  expect_equal(out$finalmod, "common-trend")
  expect_equal(out$linkfunc, "identity")
  expect_equal(
    unname(out$glm$coefficients[8]),
    0.0003801383,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with age-specific model and power5 link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "power5", pD = 0)

  expect_equal(out$finalmod, "age-specific")
  expect_equal(out$linkfunc, "power5")
  expect_equal(
    unname(out$glm$coefficients[10]),
    0.2177841,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with age-specific model and log link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "log", pD = 0)

  expect_equal(out$finalmod, "age-specific")
  expect_equal(out$linkfunc, "log")
  expect_equal(
    unname(out$glm$coefficients[3]),
    -7.833619,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with age-specific model and sqrt link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "sqrt", pD = 0)

  expect_equal(out$finalmod, "age-specific")
  expect_equal(out$linkfunc, "sqrt")
  expect_equal(
    unname(out$glm$coefficients[12]),
    0.01945599,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with age-specific model and identity link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(
    cases,
    pyr,
    nagg,
    ncase,
    linkfunc = "identity",
    pD = 0
  )

  expect_equal(out$finalmod, "age-specific")
  expect_equal(out$linkfunc, "identity")
  expect_equal(
    unname(out$glm$coefficients[18]),
    0.000161599,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("hybdproj estimate works with nba-specific model and power5 link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "power5")

  expect_equal(out$finalmod, "nba-specific")
  expect_equal(out$linkfunc, "power5")
  expect_equal(unname(out$glm$coefficients[2]), 0.1826517, tolerance = 0.000001)

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with nba-specific model and log link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "log")

  expect_equal(out$finalmod, "nba-specific")
  expect_equal(out$linkfunc, "log")
  expect_equal(
    unname(out$glm$coefficients[15]),
    -7.259301,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with nba-specific model and sqrt link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "sqrt")

  expect_equal(out$finalmod, "nba-specific")
  expect_equal(out$linkfunc, "sqrt")
  expect_equal(
    unname(out$glm$coefficients[5]),
    0.01332752,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("hybdproj estimate works with nba-specific model and identity link func", {
  cases <- matrix(2 * (1:95) + 10, nrow = 19, ncol = 5)
  pyr <- matrix(50001:50114, nrow = 19, ncol = 6)
  nagg <- 4
  ncase <- 5

  out <- hybdproj_estimate(cases, pyr, nagg, ncase, linkfunc = "identity")

  expect_equal(out$finalmod, "nba-specific")
  expect_equal(out$linkfunc, "identity")
  expect_equal(
    unname(out$glm$coefficients[18]),
    0.000161599,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})
