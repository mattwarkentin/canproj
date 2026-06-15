test_that("hybdproj prediction runs (nagg 1, average model, and power5 link func)", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "power5",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$cuttrd, 0.05)
  expect_equal(out$shortp, 0.02)
  expect_equal(out$cuttrend, 0.02)
  expect_equal(out$nopred, 1)
  expect_equal(out$noperiod, 5)
  expect_equal(out$lastperiod, 3)
  expect_equal(out$noobsper, 5)
  expect_equal(out$nototper, 6)
  expect_equal(out$noyearagg, 1)
  expect_equal(out$nocaseagp, 5)
  expect_equal(out$agrpave, integer(0))
  expect_equal(out$agrpmod, 1:19)
  expect_equal(out$linkfunc, "power5")
  expect_equal(out$projbase, 3)
  expect_equal(out$finalmod, "average")
  expect_equal(out$gofpvalue, NULL)
  expect_equal(out$glm, glm)

  expect_equal(
    out$predictions["12", "Periode 6"],
    1.860438,
    tolerance = 0.000001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.000001)
})

test_that("hybdproj prediction works with age and nba specific models w/ power5 link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.1, by = 0.005)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "power5",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "age-specific",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "power5")
  expect_equal(out$finalmod, "age-specific")

  expect_equal(
    out$predictions["7", "Periode 6"],
    39.23385,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with age and nba specific models w/ log link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.1, by = 0.005)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "log",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "nba-specific",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "log")
  expect_equal(out$finalmod, "nba-specific")

  expect_equal(
    out$predictions["17", "Periode 6"],
    85837.89,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with age and nba specific models w/ sqrt link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.1, by = 0.005)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "sqrt",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "nba-specific",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$finalmod, "nba-specific")

  expect_equal(
    out$predictions["10", "Periode 6"],
    5420.1189,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with age and nba specific models w/ identity link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      seq(from = 0.01, to = 0.1, by = 0.005)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "identity",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "age-specific",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$finalmod, "age-specific")

  expect_equal(
    out$predictions["3", "Periode 6"],
    5991.721,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with common-trend model and power5 link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      0.001
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "power5",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "common-trend",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "power5")
  expect_equal(out$finalmod, "common-trend")

  expect_equal(
    out$predictions["18", "Periode 6"],
    13.7637,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with common-trend model and log link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      0.001
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "log",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "common-trend",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "log")
  expect_equal(out$finalmod, "common-trend")

  expect_equal(
    out$predictions["11", "Periode 6"],
    56719.65,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with common-trend model and sqrt link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      0.001
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "sqrt",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "common-trend",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$finalmod, "common-trend")

  expect_equal(
    out$predictions["2", "Periode 6"],
    57.84402,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with common-trend model and identity link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01),
      0.001
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(floor(1:95 / 19) + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "identity",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "common-trend",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$finalmod, "common-trend")

  expect_equal(
    out$predictions["15", "Periode 6"],
    8217.038,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with average model and power5 link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "power5",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "power5")
  expect_equal(out$finalmod, "average")

  expect_equal(
    out$predictions["19", "Periode 6"],
    16.03648,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with average model and log link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "log",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "log")
  expect_equal(out$finalmod, "average")

  expect_equal(
    out$predictions["5", "Periode 6"],
    53198.01,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with average model and sqrt link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "sqrt",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "sqrt")
  expect_equal(out$finalmod, "average")

  expect_equal(
    out$predictions["9", "Periode 6"],
    501.0400,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with average model and identity link func", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.20, by = 0.01)
    ),
    method = "fake"
  )
  estimate_obj <- list(
    cases = matrix(2 * 1:19 + 20, nrow = 19, ncol = 5),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
    glm = glm,
    linkfunc = "identity",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:19,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(out$linkfunc, "identity")
  expect_equal(out$finalmod, "average")

  expect_equal(
    out$predictions["15", "Periode 6"],
    8017.60,
    tolerance = 0.00001
  )
  expect_snapshot_value(out$predictions, style = "json2", tolerance = 0.00001)
})

test_that("hybdproj prediction works with varying age groups", {
  glm <- list(
    coefficients = c(
      seq(from = 0.02, to = 0.11, by = 0.01)
    ),
    method = "fake"
  )

  estimate_obj <- list(
    cases = matrix(2 * 1:10 + 20, nrow = 10, ncol = 5),
    pyr = matrix(50001:50060, nrow = 10, ncol = 6),
    glm = glm,
    linkfunc = "power5",
    noyearagg = 1,
    nocaseagp = 5,
    projbase = 3,
    agrpmod = 1:10,
    agrpave = integer(0),
    finalmod = "average",
    noperiod = 5,
    lastper = 3
  )
  class(estimate_obj) <- "hybdproj.estimate"

  out <- hybdproj.prediction(estimate_obj, cuttrd = 0.05, shortp = 0.02)

  expect_equal(nrow(out$predictions), 10)
})
