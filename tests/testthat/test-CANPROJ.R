test_that("canproj can choose adpcproj during model selection", {
  cdat <- data.frame(matrix(500 + floor(1:285 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100000:100379, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop)

  expect_equal(out$method, "ADPC")
  expect_equal(out$annproj["2003", "asr"], 516.2776, tolerance = 0.00001)
  expect_equal(out$agsproj["7", "2004"], 517.1903, tolerance = 0.0001)
})

test_that("canproj can choose acproj during model selection", {
  cdat <- data.frame((row(matrix(1, 19, 15)) - col(matrix(1, 19, 15)) + 50))
  pdat <- data.frame(matrix(50000, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, pGOF = 1.0e-16, pD = 1.0e-22)

  expect_equal(out$method, "AC")
  expect_equal(out$annproj["2000", "asr"], 91.81886, tolerance = 0.00001)
  expect_equal(out$agsproj["9", "2002"], 83.03031, tolerance = 0.0001)
})

test_that("canproj can choose hybdproj during model selection", {
  cdat <- data.frame(matrix(
    c(rep(1, 19), rep(5, 19), rep(2, 19)),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop)

  expect_equal(out$method, "Hybrid")
  expect_equal(out$annproj["2001", "asr"], 26.076498, tolerance = 0.00001)
  expect_equal(out$agsproj["3", "2004"], 26.09135, tolerance = 0.0001)
})

test_that("canproj can choose ave5proj during model selection", {
  cdat <- data.frame(matrix(
    rep(c(rep(0, 52), rep(1, 5)), 5),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop)

  expect_equal(out$method, "average5")

  expect_equal(out$annproj["2003", "asr"], 0.492162, tolerance = 0.00001)
  expect_snapshot_value(out$annproj, style = "json2", tolerance = 0.00001)

  expect_equal(out$agsproj["15", "2002"], 3.905487, tolerance = 0.00001)
  expect_snapshot_value(out$agsproj, style = "json2", tolerance = 0.00001)
})

test_that("canproj can specify nordpred method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "nordpred")

  expect_equal(out$method, "nordpred")
  expect_equal(out$annproj["2003", "asr"], 318.0060, tolerance = 0.00001)
  expect_equal(out$agsproj["14", "2004"], 318.9471, tolerance = 0.0001)
})

test_that("canproj can specify adpc-nb method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "adpc-nb")

  expect_equal(out$method, "adpc-nb")
  expect_equal(out$annproj["2001", "asr"], 316.0908, tolerance = 0.00001)
  expect_equal(out$agsproj["17", "2001"], 316.0651, tolerance = 0.0001)
})

test_that("canproj can specify ac-poi method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "ac-poi")

  expect_equal(out$method, "ac-poi")
  expect_equal(out$annproj["2002", "asr"], 312.8873, tolerance = 0.00001)
  expect_equal(out$agsproj["9", "2003"], 318.0188, tolerance = 0.00001)
})

test_that("canproj can specify ac-nb method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "ac-nb")

  expect_equal(out$method, "ac-nb")
  expect_equal(out$annproj["2003", "asr"], 313.0127, tolerance = 0.00001)
  expect_equal(out$agsproj["18", "2000"], 315.1104, tolerance = 0.00001)
})

test_that("canproj can specify age-trd-nb method", {
  set.seed(54638)
  cdat <- data.frame(matrix(
    rnbinom(19 * 15, size = 10, mu = 320),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "age-trd-nb")

  expect_equal(out$method, "a-s-nb")
  expect_equal(out$annproj["2004", "asr"], 380.3885, tolerance = 0.00001)
  expect_equal(out$agsproj["15", "2001"], 285.9604, tolerance = 0.00001)
})

test_that("canproj can specify age-trd-poi method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "age-trd-poi")

  expect_equal(out$method, "a-s-poi")
  expect_equal(out$annproj["2002", "asr"], 317.0568, tolerance = 0.00001)
  expect_equal(out$agsproj["3", "2003"], 318.0399, tolerance = 0.00001)
})

test_that("canproj can specify com-trd method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "com-trd")

  expect_equal(out$method, "c-t")
  expect_equal(out$annproj["2001", "asr"], 316.0946, tolerance = 0.00001)
  expect_equal(out$agsproj["7", "2000"], 315.1404, tolerance = 0.00001)
})

test_that("canproj can specify age-only method", {
  cdat <- data.frame(matrix(300 + floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "age-only")

  expect_equal(out$method, "average")
  expect_equal(out$annproj["2004", "asr"], 307.5638, tolerance = 0.00001)
  expect_equal(out$agsproj["2", "2002"], 307.5848, tolerance = 0.0001)
})

test_that("canproj can specify ave5 method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  standpop <- stdpop_Canada_2021

  out <- canproj(cdat, pdat, 2000, standpop, methods = "ave5")

  expect_equal(out$method, "average5")
  expect_equal(out$annproj["2001", "asr"], 244.13305, tolerance = 0.00001)
  expect_equal(out$agsproj["14", "2003"], 311.2468, tolerance = 0.0001)
})

test_that("canproj works with varying age groups", {
  cdat <- data.frame(matrix(95, nrow = 10, ncol = 15))
  pdat <- data.frame(matrix(10000:10199, nrow = 10, ncol = 20))
  standpop <- StandardPopulation("dummy", rep("a", 10), rep(0.1, 10))
  out <- canproj(cdat, pdat, 2000, standpop)

  expect_equal(nrow(out$agsproj), 10)
  expect_equal(nrow(out$out$predictions), 10)
})
