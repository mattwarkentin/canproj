test_that("Canada 2021 standard population works", {
  expect_equal(
    stdpop_Canada_2021@weights,
    c(
      0.049762,
      0.054199,
      0.055283,
      0.053870,
      0.062877,
      0.069944,
      0.070712,
      0.069120,
      0.065514,
      0.062285,
      0.063488,
      0.070359,
      0.068330,
      0.058224,
      0.048268,
      0.033273,
      0.022004,
      0.013635,
      0.008853
    )
  )

  expect_equal(
    stdpop_Canada_2021@strata,
    c(
      "0 to 4 years",
      "5 to 9 years",
      "10 to 14 years",
      "15 to 19 years",
      "20 to 24 years",
      "25 to 29 years",
      "30 to 34 years",
      "35 to 39 years",
      "40 to 44 years",
      "45 to 49 years",
      "50 to 54 years",
      "55 to 59 years",
      "60 to 64 years",
      "65 to 69 years",
      "70 to 74 years",
      "75 to 79 years",
      "80 to 84 years",
      "85 to 89 years",
      "90 years and over"
    )
  )
})

test_that("World standard population works", {
  expect_equal(
    stdpop_WHO_2000_2025@weights,
    c(
      0.08857,
      0.08687,
      0.08597,
      0.08467,
      0.08217,
      0.07927,
      0.07607,
      0.07148,
      0.06588,
      0.06038,
      0.05368,
      0.04548,
      0.03719,
      0.02959,
      0.02209,
      0.01519,
      0.00910,
      0.00440,
      0.00195
    )
  )
})

test_that("Canada 2011 standard population works", {
  expect_equal(
    stdpop_Canada_2011@weights,
    c(
      0.055325,
      0.052736,
      0.055857,
      0.065121,
      0.068506,
      0.068967,
      0.067743,
      0.066176,
      0.069477,
      0.079209,
      0.078366,
      0.068519,
      0.059696,
      0.044613,
      0.033573,
      0.026761,
      0.020406,
      0.012420,
      0.006529
    )
  )
})

test_that("standard population class catches invalid inputs", {
  expect_error(
    StandardPopulation(name = "Canada 2021", strata = "all ages"),
    "@weights is required."
  )

  expect_error(
    StandardPopulation(name = "Canada 2021", weights = c(0.5, 0.5)),
    "@strata is required."
  )

  expect_error(
    StandardPopulation(weights = c(0.5, 0.5), strata = "all ages"),
    "@name is required."
  )

  expect_error(
    StandardPopulation(
      name = c("Canada 2021", "other"),
      strata = "all ages",
      weights = c(0.5, 0.5)
    ),
    "`name` must have length 1."
  )

  expect_error(
    StandardPopulation(
      name = c("Canada 2021"),
      strata = c(1:3),
      weights = c(0.5, 0.5, 0.5)
    ),
    "Must sum to 1."
  )

  expect_error(
    StandardPopulation(
      name = c("Canada 2021"),
      strata = "all ages",
      weights = c(0.5, 0.5)
    ),
    "`strata` and `weights` must be the same length."
  )
})

test_that("standard popultation class can be instantiated", {
  pop <- StandardPopulation(
    "Can",
    c("0-29", "30-59", "60+"),
    c(0.35, 0.35, 0.3),
    list("dummy population")
  )

  expect_equal(pop@name, "Can")
  expect_equal(pop@strata, c("0-29", "30-59", "60+"))
  expect_equal(pop@weights, c(0.35, 0.35, 0.3))
  expect_equal(pop@metadata, list("dummy population"))
})

test_that("as.data.frame behaves correctly for standard population class", {
  pop <- StandardPopulation(
    "Can",
    c("0-29", "30-59", "60+"),
    c(0.35, 0.35, 0.3),
    list("dummy population")
  )

  expected <- data.frame(
    strata = c("0-29", "30-59", "60+"),
    weights = c(0.35, 0.35, 0.30)
  )

  expect_equal(as.data.frame(pop), expected)
})
