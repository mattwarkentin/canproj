test_that("acproj works", {
  cdat <- matrix(floor((0:284) / 19) + 5, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  output <- acproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    startage = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    pGOF = 0.05,
    linkfunc = "power5"
  )

  expect_s3_class(output, "acproj")
  expect_equal(output$nopred, 1)
  expect_equal(output$noperiod, 3)
  expect_equal(output$cuttrd, 0.05)
  expect_equal(output$shortp, 0.01)
  expect_equal(output$cuttrend, 0.01)
  expect_equal(output$distribution, "Poisson")
  expect_equal(output$startage, 1)
  expect_s3_class(output$glm, "glm")

  expect_equal(output$gofpvalue, 0.9997187, tolerance = 0.0000001)
  expect_equal(output$predictions["4", "4"], 127.04251, tolerance = 0.00001)
  expect_equal(output$pyr["10", "3"], 51185)
  expect_equal(
    unname(output$glm$coefficients[14]),
    0.3444149,
    tolerance = 0.000001
  )

  expect_snapshot_value(
    output$predictions,
    style = "json2",
    tolerance = 0.000001
  )
  expect_snapshot_value(output$pyr, style = "json2", tolerance = 0.000001)
  expect_snapshot_value(
    output$glm$coefficients,
    style = "json2",
    tolerance = 0.000001
  )
})

test_that("acproj must have at least 15 years in cases", {
  cdat <- matrix(1, nrow = 19, ncol = 10)
  pdat <- matrix(1000, nrow = 19, ncol = 15)

  expect_error(
    acproj(cdat, pdat),
    "Minimum number of period is 3 (15 years) in \"cases\"",
    fixed = TRUE
  )
})

test_that("acproj moves start stage according to n5case when NULL", {
  cdat <- matrix(rep(1:19, 15), nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  output <- acproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 6,
    startage = NULL,
    cuttrd = 0.04,
    shortp = 0,
    pGOF = 0.05,
    linkfunc = "power5"
  )

  expect_equal(output$startage, 2)
})

test_that("acproj works with varying age groups", {
  cdat <- matrix(floor((0:149) / 19) + 5, nrow = 10, ncol = 15)
  pdat <- matrix(10000:10199, nrow = 10, ncol = 20)

  output <- acproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    startage = 1,
    cuttrd = 0.05,
    shortp = 0.01,
    pGOF = 0.05,
    linkfunc = "power5"
  )

  expect_equal(nrow(output$predictions), 10)
  expect_equal(output$predictions["4", "4"], 73.09689, tolerance = 0.00001)
})

test_that("acproj getproj works with a standard population", {
  cdat <- matrix(floor((0:284) / 19) + 5, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(acproj_obj) <- "acproj"

  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)
  )

  out <- acproj_get_projections(cdat, pdat, 2000, acproj_obj, standpop = stdpop)

  expect_equal(out["2002", "asr"], 217.72721, tolerance = 0.00001)
  expect_equal(out["1990", "case"], 190, tolerance = 0.00001)

  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("acproj getproj works without a standard population", {
  cdat <- matrix(floor((0:284) / 19) + 5, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  acproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(acproj_obj) <- "acproj"

  out <- acproj_get_projections(cdat, pdat, 2000, acproj_obj)

  expect_equal(out["5", "1988"], 79.51496, tolerance = 0.00001)
  expect_equal(out["12", "2003"], 227.5242, tolerance = 0.00001)

  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("acproj plot works", {
  skip_if_not_installed("vdiffr")

  cdat <- matrix(floor((0:284) / 19) + 5, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)
  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)
  )

  output <- acproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    startage = 1,
    cuttrd = 0.04,
    shortp = 0,
    pGOF = 0.05,
    linkfunc = "power5"
  )

  vdiffr::expect_doppelganger(
    "acproj plot",
    function() plot(output, cdat, pdat, 2000, stdpop)
  )
})
