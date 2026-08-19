test_that("hybdproj works", {
  cdat <- data.frame(matrix(
    c(rep(1, 19), rep(5, 19), rep(2, 19)),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- hybdproj(cdat, pdat, stdpop_Canada_2021)

  expect_equal(out$cuttrd, 0.04)
  expect_equal(out$shortp, 0)
  expect_equal(out$cuttrend, rep(0, 5))
  expect_equal(out$nopred, 5)
  expect_equal(out$noperiod, 6)
  expect_equal(out$lastperiod, 6)
  expect_equal(out$noobsper, 15)
  expect_equal(out$nototper, 20)
  expect_equal(out$noyearagg, 1)
  expect_equal(out$nocaseagp, 1)
  expect_equal(out$agrpave, integer(0))
  expect_equal(out$agrpmod, 1:19)
  expect_equal(out$linkfunc, "power5")
  expect_equal(out$projbase, 6)
  expect_equal(out$finalmod, "average")
  expect_equal(out$gofpvalue, 0.03994028, tolerance = 0.000001)

  expect_equal(out$predictions["3", "18"], 2.693932, tolerance = 0.000001)
  expect_s3_class(out$glm, "glm")
})

test_that("hybdproj works with varying age groups", {
  cdat <- data.frame(matrix(
    c(rep(1, 10), rep(5, 10), rep(2, 10)),
    nrow = 10,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10199, nrow = 10, ncol = 20))
  pop <- StandardPopulation("dummy", rep("a", 10), rep(0.1, 10))

  out <- hybdproj(cdat, pdat, pop)

  expect_equal(length(out$glm$coefficients), 10)
  expect_equal(out$agrpmod, 1:10)
  expect_equal(nrow(out$predictions), 10)
})

test_that("hybdproj plotting works", {
  cdat <- data.frame((row(matrix(1, 19, 15)) - col(matrix(1, 19, 15)) + 50))
  pdat <- data.frame(matrix(50000, nrow = 19, ncol = 20))

  out <- hybdproj(cdat, pdat, stdpop_Canada_2021)

  vdiffr::expect_doppelganger("acproj plot", function() {
    plot(out, cdat, pdat, 2020, stdpop_Canada_2021)
  })
})
