test_that("acproj getpred requires acproj.object class", {
  object <- list(1:6)

  expect_error(
    acproj.getpred(object, TRUE),
    "Variable \"acproj.object\" must be of type \"acproj\"",
    fixed = TRUE
  )
})

test_that("acproj getpred requires incidence if using a standard population", {
  acproj_obj <- list(1:5)
  class(acproj_obj) <- "acproj"

  stdpop <- 1:19

  expect_error(
    acproj.getpred(acproj_obj, incidence = F, standpop = stdpop),
    "\"standpop\" should only be used with incidence predictions (incidence=T)",
    fixed = TRUE
  )
})

test_that("acproj getpred requires 'by age' = T when using standard population", {
  acproj_obj <- list(1:5)
  class(acproj_obj) <- "acproj"

  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expect_error(
    acproj.getpred(acproj_obj, standpop = stdpop, byage = TRUE),
    "\"standpop\" is only valid for \"byage=F\"",
    fixed = TRUE
  )
})

test_that("acproj getpred can specify age groups to use", {
  acproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(acproj_obj) <- "acproj"

  expected <- data.frame(
    X1 = c(49.99500049995, 898.382910760631),
    X2 = c(998.003992015968, 1843.18023313739),
    X3 = c(1942.42454427732, 2784.40731901352),
    X4 = c(2883.27699343806, 3722.08436724566)
  )
  attr(expected, "row.names") <- as.integer(c(2, 19))

  expect_equal(
    acproj.getpred(acproj_obj, agegroups = c(2, 19)),
    expected,
    tolerance = 0.0001
  )
})

test_that("acproj getpred works with a standard population", {
  acproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(acproj_obj) <- "acproj"

  stdpop <- StandardPopulation(
    "a",
    rep("a", 19),
    c(rep(0.04, 10), rep(0.07, 8), 0.04)
  )

  expected <- c(
    X1 = 503.354104249281,
    X2 = 1449.64643772789,
    X3 = 2392.36006379238,
    X4 = 3331.51524517502
  )

  expect_equal(
    acproj.getpred(acproj_obj, standpop = stdpop, byage = FALSE),
    expected,
    tolerance = 0.0001
  )
})

test_that("acproj getpred works when excluding observed values", {
  acproj_obj <- list(
    predictions = matrix(0:75 * 25, nrow = 19, ncol = 4),
    pyr = matrix(seq(from = 50000, to = 50375, by = 5), nrow = 19, ncol = 4),
    nopred = 1
  )
  class(acproj_obj) <- "acproj"

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
    acproj.getpred(acproj_obj, excludeobs = TRUE),
    expected,
    tolerance = 0.0001
  )
})
