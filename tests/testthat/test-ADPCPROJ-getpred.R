test_that("adpcproj getpred requires object class adpcproj", {
  object <- list(
    predictions = matrix(5, 19, 4),
    pyr = matrix(10, 19, 4)
  )

  expect_error(
    adpcproj_get_predictions(object),
    "Variable \"adpcproj.object\" must be of type \"adpcproj\"",
    fixed = TRUE
  )
})

test_that("adpcproj getpred standpop requires incidence prediction", {
  adpc_obj <- list(
    predictions = matrix(5, 19, 4),
    pyr = matrix(10, 19, 4)
  )
  class(adpc_obj) <- "adpcproj"

  stdpop <- 1:19

  expect_error(
    adpcproj_get_predictions(adpc_obj, incidence = FALSE, standpop = stdpop),
    "\"standpop\" should only be used with incidence predictions (incidence=T)",
    fixed = TRUE
  )
})

test_that("adpcproj getpred standpop must be same length as age groups", {
  adpc_obj <- list(
    predictions = matrix(5, 19, 4),
    pyr = matrix(10, 19, 4)
  )
  class(adpc_obj) <- "adpcproj"

  stdpop <- StandardPopulation("a", rep("a", 2), c(0.5, 0.5))

  expect_error(
    adpcproj_get_predictions(adpc_obj, agegroups = 4:9, standpop = stdpop),
    "\"standpop\" must be the same length as \"agegroups\"",
    fixed = TRUE
  )
})

test_that("adpcproj getpred standpop requires 'by age' = F", {
  adpc_obj <- list(
    predictions = matrix(5, 19, 4),
    pyr = matrix(10, 19, 4)
  )
  class(adpc_obj) <- "adpcproj"

  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expect_error(
    adpcproj_get_predictions(adpc_obj, byage = TRUE, standpop = stdpop),
    "\"standpop\" is only valid for \"byage=F\"",
    fixed = TRUE
  )
})

test_that("adpcproj getpred can select specific age groups", {
  adpcproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(adpcproj_obj) <- "adpcproj"

  expected <- data.frame(
    X1 = c(99.98000, 349.75517),
    X2 = c(1047.7996, 1296.6288),
    X3 = c(1992.0319, 2239.9204),
    X4 = c(2932.6971, 3179.6502)
  )
  attr(expected, "row.names") <- as.integer(c(3, 8))

  expect_equal(
    adpcproj_get_predictions(adpcproj_obj, agegroups = c(3, 8)),
    expected,
    tolerance = 0.0001
  )
})

test_that("adpcproj getpred works with a standard population", {
  adpcproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(adpcproj_obj) <- "adpcproj"

  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expected <- c(X1 = 503.3541, X2 = 1449.6464, X3 = 2392.3601, X4 = 3331.5152)

  expect_equal(
    adpcproj_get_predictions(adpcproj_obj, standpop = stdpop),
    expected,
    tolerance = 0.0001
  )
})

test_that("adpcproj getpred works when excluding observed values", {
  adpcproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(adpcproj_obj) <- "adpcproj"

  expected <- data.frame(
    X4 = c(
      2833.84707169136,
      2883.27699343806,
      2932.69708718561,
      2982.10735586481,
      3031.50780240533,
      3080.89842973564,
      3130.27924078307,
      3179.65023847377,
      3229.01142573274,
      3278.36280548381,
      3327.70438064965,
      3377.03615415177,
      3426.35812891052,
      3475.67030784508,
      3524.9726938735,
      3574.26528991263,
      3623.54809887819,
      3672.82112368473,
      3722.08436724566
    )
  )
  rownames(expected) <- 1:19

  expect_equal(
    adpcproj_get_predictions(adpcproj_obj, excludeobs = TRUE),
    expected,
    tolerance = 0.0000001
  )
})
