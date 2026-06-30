test_that("hybdproj getpred allows for selecting age groups", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
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

  expected <- data.frame(
    X1 = c(69.9972001119955, 69.9860027994401),
    X2 = c(119.949621159113, 119.9304403446),
    X3 = c(169.86410871303, 169.836956521739),
    X4 = c(219.740705966959, 219.705594503365),
    X5 = c(69.8909700866648, 69.8798067324202),
    X6 = c(119.76765075753, 119.748528091009)
  )
  attr(expected, "row.names") <- as.integer(c(2, 10))

  expect_equal(
    hybdproj_get_predictions(hybd_obj, incidence = T, agegroups = c(2, 10)),
    expected,
    tolerance = 0.000001
  )
})

test_that("hybdproj getpred works with a standard population", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
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

  stdpop <- StandardPopulation(
    "dummy",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expected <- c(
    X1 = 69.9844922133113,
    X2 = 119.927852734417,
    X3 = 169.833293522694,
    X4 = 219.700857747862,
    X5 = 69.878300725886,
    X6 = 119.745948322627
  )

  expect_equal(
    hybdproj_get_predictions(hybd_obj, standpop = stdpop),
    expected,
    tolerance = 0.000001
  )
})

test_that("hybdproj getpred can track by case numbers", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
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

  expected <- matrix(floor(0:75 / 19) * 25 + 35, 19, 4)

  expect_equal(
    hybdproj_get_predictions(hybd_obj, incidence = F),
    expected,
    tolerance = 0.000001
  )
})

test_that("hybdproj getpred works without observed values", {
  hybd_obj <- list(
    predictions = matrix(floor(0:75 / 19) * 25 + 35, nrow = 19, ncol = 4),
    pyr = matrix(50001:50114, nrow = 19, ncol = 6),
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

  expected <- data.frame(
    X6 = c(
      119.770041520281,
      119.76765075753,
      119.765260090223,
      119.762869518354,
      119.760479041916,
      119.758088660905,
      119.755698375314,
      119.753308185139,
      119.750918090372,
      119.748528091009,
      119.746138187043,
      119.74374837847,
      119.741358665283,
      119.738969047476,
      119.736579525045,
      119.734190097982,
      119.731800766284,
      119.729411529942,
      119.727022388953
    )
  )
  rownames(expected) <- 1:19

  expect_equal(
    hybdproj_get_predictions(hybd_obj, excludeobs = T),
    expected,
    tolerance = 0.000001
  )
})
