test_that("canproj can choose adpcproj during model selection", {
  cdat <- data.frame(matrix(floor(0:284 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000)

  expect_equal(out$method, "ADPC")
  expect_equal(out$annproj["2003", "asr"], 203.65801, tolerance = 0.00001)
  expect_equal(out$agsproj["30-34", "2004"], 218.1044, tolerance = 0.0001)
})

test_that("canproj can choose acproj during model selection", {
  cdat <- data.frame((row(matrix(1, 19, 15)) - col(matrix(1, 19, 15)) + 50))
  pdat <- data.frame(matrix(50000, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000, pGOF = 1.0e-16, pD = 1.0e-22)

  expect_equal(out$method, "AC")
  expect_equal(out$annproj["2000", "asr"], 91.81886, tolerance = 0.00001)
  expect_equal(out$agsproj["40-44", "2002"], 83.03031, tolerance = 0.0001)
})

test_that("canproj can choose hybdproj during model selection", {
  cdat <- data.frame(matrix(
    c(rep(1, 19), rep(5, 19), rep(2, 19)),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000)

  expect_equal(out$method, "Hybrid")
  expect_equal(out$annproj["2001", "asr"], 26.076498, tolerance = 0.00001)
  expect_equal(out$agsproj["10-14", "2004"], 26.09135, tolerance = 0.0001)
})

test_that("canproj can choose ave5proj during model selection", {
  cdat <- data.frame(matrix(
    rep(c(rep(0, 52), rep(1, 5)), 5),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000)

  expect_equal(out$method, "average5")

  expect_equal(out$annproj["2003", "asr"], 0.492162, tolerance = 0.00001)
  expect_snapshot_value(out$annproj, style = "json2", tolerance = 0.00001)

  expect_equal(out$agsproj["70-74", "2002"], 3.905487, tolerance = 0.00001)
  expect_snapshot_value(out$agsproj, style = "json2", tolerance = 0.00001)
})

test_that("canproj can specify nordpred method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000, methods = "nordpred")

  expect_equal(out$method, "nordpred")
  expect_equal(out$annproj["2003", "asr"], 404.35309, tolerance = 0.00001)
  expect_equal(out$agsproj["65-69", "2004"], 553.8477, tolerance = 0.0001)
})

test_that("canproj can specify adpc-nb method", {
  cdat <- data.frame(matrix(1:285, nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000, methods = "adpc-nb")

  expect_equal(out$method, "adpc-nb")
  expect_equal(out$annproj["2001", "asr"], 47.25288, tolerance = 0.00001)
  expect_equal(out$agsproj["80-84", "2001"], 48.88024, tolerance = 0.0001)
})

test_that("canproj can specify ac-poi method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000, methods = "ac-poi")

  expect_equal(out$method, "ac-poi")
  expect_equal(out$annproj["2002", "asr"], 399.25027, tolerance = 0.00001)
  expect_equal(out$agsproj["40-44", "2003"], 488.12615, tolerance = 0.00001)
})

test_that("canproj can specify ac-nb method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))

  out <- canproj(cdat, pdat, 2000, methods = "ac-nb")

  expect_equal(out$method, "ac-nb")
  expect_equal(out$annproj["2003", "asr"], 466.08970, tolerance = 0.00001)
  expect_equal(out$agsproj["85-89", "2000"], 570.75521, tolerance = 0.00001)
})

test_that("canproj can specify age-trd-nb method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  out <- canproj(cdat, pdat, 2000, methods = "age-trd-nb")

  expect_equal(out$method, "a-s-nb")
  expect_equal(out$annproj["2004", "asr"], 509.25829, tolerance = 0.00001)
  expect_equal(out$agsproj["70-74", "2001"], 505.1084, tolerance = 0.00001)
})

test_that("canproj can specify age-trd-poi method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  out <- canproj(cdat, pdat, 2000, methods = "age-trd-poi")

  expect_equal(out$method, "a-s-poi")
  expect_equal(out$annproj["2002", "asr"], 410.22099, tolerance = 0.00001)
  expect_equal(out$agsproj["10-14", "2003"], 311.9490, tolerance = 0.00001)
})

test_that("canproj can specify com-trd method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  out <- canproj(cdat, pdat, 2000, methods = "com-trd")

  expect_equal(out$method, "c-t")
  expect_equal(out$annproj["2001", "asr"], 374.76110, tolerance = 0.00001)
  expect_equal(out$agsproj["30-34", "2000"], 310.6442, tolerance = 0.00001)
})

test_that("canproj can specify age-only method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  out <- canproj(cdat, pdat, 2000, methods = "age-only")

  expect_equal(out$method, "average")
  expect_equal(out$annproj["2004", "asr"], 501.74981, tolerance = 0.00001)
  expect_equal(out$agsproj["5-9", "2002"], 284.9060, tolerance = 0.0001)
})

test_that("canproj can specify ave5 method", {
  cdat <- data.frame(matrix(11:29 * floor(19:303 / 19), nrow = 19, ncol = 15))
  pdat <- data.frame(matrix(100001:100380, nrow = 19, ncol = 20))
  out <- canproj(cdat, pdat, 2000, methods = "ave5")

  expect_equal(out$method, "average5")
  expect_equal(out$annproj["2001", "asr"], 244.13305, tolerance = 0.00001)
  expect_equal(out$agsproj["65-69", "2003"], 311.2468, tolerance = 0.0001)
})

test_that("canproj plot works", {
  cdat <- matrix(floor((0:284) / 19) + 5, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)
  stdpop <- c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)

  out <- canproj(cdat, pdat, 2000)

  vdiffr::expect_doppelganger("canproj plot", function() plot(out, stdpop))
})
