test_that("adpcproj works with all variables", {
  cdat <- matrix(floor(0:284 / 19) + 20, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  out <- adpcproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    noperiods = 3,
    recent = TRUE,
    startestage = 1,
    newcohort = NULL,
    pGOF = 0.1,
    cuttrd = 0.03,
    shortp = 0.02,
    linkfunc = "power5"
  )

  expect_s3_class(out, "adpcproj")
  expect_equal(out$nopred, 1)
  expect_equal(out$noperiod, 3)
  expect_equal(out$gofpvalue, 1)
  expect_equal(out$recent, TRUE)
  expect_equal(out$pvaluerecent, 0.4815863, tolerance = 0.000001)
  expect_equal(out$cuttrd, 0.03)
  expect_equal(out$shortp, 0.02)
  expect_equal(out$cuttrend, 0.02)
  expect_equal(out$distribution, "Poisson")
  expect_equal(out$startuseage, 1)
  expect_equal(out$startestage, 1)

  expect_equal(out$predictions["5-9", "4"], 188.0918, tolerance = 0.0001)
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)

  expect_equal(out$pyr["12", "4"], 51670)
  expect_snapshot_value(out$pyr, style = "json2", tolerance = 0.0001)

  expect_equal(unname(out$glm$coefficients[7]), 0.2830023, tolerance = 0.000001)
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("adpcproj assigns start stage when NULL", {
  cdat <- data.frame(matrix(0:18 * 0:284, nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- adpcproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    noperiods = 3,
    recent = FALSE,
    startestage = NULL,
    newcohort = NULL,
    pGOF = 0.1,
    cuttrd = 0.03,
    shortp = 0.02,
    linkfunc = "power5"
  )

  expect_equal(out$startuseage, 2)
  expect_equal(out$startestage, 2)

  expect_length(out$glm$coefficients, 40)
  expect_equal(
    unname(out$glm$coefficients[9]),
    0.3612491,
    tolerance = 0.000001
  )
  expect_snapshot_value(
    out$glm$coefficients,
    style = "json2",
    tolerance = 0.00001
  )
})

test_that("adpcproj requires at least 15 years in cases", {
  cdat <- data.frame(matrix(0:18 * 0:189, nrow = 19, ncol = 10))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  expect_error(
    adpcproj(cdat, pdat),
    "Minimum number of period is 3 (15 years) in \"cases\"",
    fixed = TRUE
  )
})

test_that("adpcproj determines number of periods when NULL", {
  cdat <- matrix(floor(0:284 / 19) + 20, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  out <- adpcproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    noperiods = NULL,
    recent = FALSE,
    startestage = 1,
    newcohort = NULL,
    pGOF = 0.1,
    cuttrd = 0.03,
    shortp = 0.02,
    linkfunc = "power5"
  )

  expect_equal(out$noperiod, 3)
})

test_that("adpcproj decides between recent or whole trend when NULL", {
  cdat <- matrix(floor(0:284 / 19) + 20, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)

  out <- adpcproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    noperiods = NULL,
    recent = NULL,
    startestage = 1,
    newcohort = NULL,
    pGOF = 0.1,
    cuttrd = 0.03,
    shortp = 0.02,
    linkfunc = "power5"
  )

  expect_equal(out$recent, FALSE)
})

test_that("adpcproj getproj works with standard population", {
  cdat <- matrix(floor((0:284) / 19), nrow = 19, ncol = 15)
  pdat <- matrix(50000:50379, nrow = 19, ncol = 20)

  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(adpcproj_obj) <- "adpcproj"

  stdpop <- c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)

  out <- adpcproj.getproj(cdat, pdat, 2000, adpcproj_obj, standpop = stdpop)

  expect_equal(out["1994", "asr"], 17.935293, tolerance = 0.00001)
  expect_equal(out["2000", "case"], 1893, tolerance = 0.00001)

  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("adpcproj getproj works without standard population", {
  cdat <- matrix(floor((0:284) / 19), nrow = 19, ncol = 15)
  pdat <- matrix(50000:50379, nrow = 19, ncol = 20)

  adpcproj_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50190, to = 50565, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(adpcproj_obj) <- "adpcproj"

  out <- adpcproj.getproj(cdat, pdat, 2000, adpcproj_obj)

  expect_equal(out["90+", "1987"], 3.995525, tolerance = 0.00001)
  expect_equal(out["25-29", "2000"], 198.1467, tolerance = 0.00001)

  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("adpcproj plot works", {
  cdat <- matrix(floor(0:284 / 19) + 20, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)
  stdpop <- c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)

  out <- adpcproj(
    cdat,
    pdat,
    projfor = "incidence",
    n5case = 5,
    noperiods = 3,
    recent = TRUE,
    startestage = 1,
    newcohort = NULL,
    pGOF = 0.1,
    cuttrd = 0.03,
    shortp = 0.02,
    linkfunc = "power5"
  )

  vdiffr::expect_doppelganger(
    "adpcproj plot",
    function() plot(out, cdat, pdat, 2000, stdpop)
  )
})
