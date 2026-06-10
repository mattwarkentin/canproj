test_that("hybdproj works", {
  cdat <- data.frame(matrix(
    c(rep(1, 19), rep(5, 19), rep(2, 19)),
    nrow = 19,
    ncol = 15
  ))
  pdat <- data.frame(matrix(10000:10379, nrow = 19, ncol = 20))

  out <- hybdproj(cdat, pdat)

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

  expect_equal(out$predictions["10-14", "18"], 2.693932, tolerance = 0.000001)
  expect_s3_class(out$glm, "glm")
})

test_that("hybdproj plot works", {
  skip_if_not_installed("vdiffr")

  cdat <- matrix(floor(0:284 / 19) + 20, nrow = 19, ncol = 15)
  pdat <- matrix(10000:10379, nrow = 19, ncol = 20)
  stdpop <- c(rep(0.05, 15), 0.06, 0.06, 0.06, 0.07)

  out <- hybdproj(
    cdat,
    pdat,
    projfor = "incidence",
    nagg = 1,
    ncase = 5,
    cuttrd = 0.05,
    shortp = 0.01,
    linkfunc = "power5",
    pD = 0.05,
    pGOF = 0.05
  )

  vdiffr::expect_doppelganger(
    "hybdproj plot",
    function() plot(out, cdat, pdat, 2000, stdpop)
  )
})
