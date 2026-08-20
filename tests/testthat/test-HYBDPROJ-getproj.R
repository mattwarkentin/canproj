test_that("hybdproj getproj works with non-average model", {
  cdat <- matrix(1:228, 19, 12)
  pdat <- matrix(50000, 19, 16)

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
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

  out <- get_projections(
    hybd_obj,
    cdat = cdat,
    pdat = pdat,
    startp = 2000
  )

  expect_equal(out["11", "2002"], 57.41385, tolerance = 0.00001)
  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj getproj works with non-average model and standpop", {
  cdat <- matrix(1:228, 19, 12)
  pdat <- matrix(50000, 19, 16)

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
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

  stdpop <- StandardPopulation(
    "dummy",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  suppressWarnings(
    out <- get_projections(
      hybd_obj,
      cdat = cdat,
      pdat = pdat,
      startp = 2000,
      standpop = stdpop
    )
  )

  expect_equal(out["2000", "asr"], 107.37800, tolerance = 0.00001)
  expect_equal(out["2000", "case"], 1020)
  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj getproj works with 5-year average method", {
  cdat <- data.frame(matrix(1:228, 19, 12))
  pdat <- data.frame(matrix(50000, 19, 16))

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "average",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  out <- get_projections(
    hybd_obj,
    cdat = cdat,
    pdat = pdat,
    startp = 2000,
    ave5 = TRUE
  )

  expect_equal(out["17", "2003"], 376)
  expect_snapshot_value(out, style = "json2")
})

test_that("hybdproj getproj works with 5-year average method and standpop", {
  cdat <- data.frame(matrix(1:228, 19, 12))
  pdat <- data.frame(matrix(50000, 19, 16))

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "average",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  stdpop <- StandardPopulation(
    "dummy",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  suppressWarnings(
    out <- get_projections(
      hybd_obj,
      cdat = cdat,
      pdat = pdat,
      startp = 2000,
      ave5 = TRUE,
      standpop = stdpop
    )
  )

  expect_equal(out["2001", "asr"], 364.16, tolerance = 0.001)
  expect_equal(out["2001", "case"], 3439)
  expect_snapshot_value(out, style = "json2", tolerance = 0.001)
})

test_that("hybdproj getproj works with average method", {
  cdat <- data.frame(matrix(1:228, 19, 12))
  pdat <- data.frame(matrix(50000, 19, 16))

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "average",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  out <- get_projections(hybd_obj, cdat = cdat, pdat = pdat, startp = 2000)

  expect_equal(out["7", "2002"], 57.41844, tolerance = 0.00001)
  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj getproj works with average method and standpop", {
  cdat <- data.frame(matrix(1:228, 19, 12))
  pdat <- data.frame(matrix(50000, 19, 16))

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
    noobsper = 3,
    nopred = 1,
    noyearagg = 4,
    noperiod = 12,
    shortp = 0.02,
    cuttrd = 0.05,
    projbase = 13,
    nocaseagp = 1,
    finalmod = "average",
    linkfunc = "power5",
    agrpmod = "1,2,3,4,5",
    agrpave = c(1, 2, 3, 6, 7, 8),
    gofpvalue = 0.812345
  )
  class(hybd_obj) <- "hybdproj"

  stdpop <- StandardPopulation(
    "dummy",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  suppressWarnings(
    out <- get_projections(
      hybd_obj,
      cdat = cdat,
      pdat = pdat,
      startp = 2000,
      standpop = stdpop
    )
  )

  expect_equal(out["2003", "asr"], 32.43164, tolerance = 0.00001)
  expect_equal(out["2003", "case"], 308)
  expect_snapshot_value(out, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj getproj triggers apc warning", {
  cdat <- matrix(1:228, 19, 12)
  pdat <- matrix(50000, 19, 16)

  hybd_obj <- list(
    predictions = matrix(floor(0:56 / 19) * 25 + 35, nrow = 19, ncol = 3),
    pyr = matrix(50001:50076, nrow = 19, ncol = 4),
    nototper = 3,
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

  stdpop <- StandardPopulation(
    "dummy",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expect_warning(
    out <- get_projections(
      hybd_obj,
      cdat = cdat,
      pdat = pdat,
      startp = 2000,
      standpop = stdpop
    )
  )
})
