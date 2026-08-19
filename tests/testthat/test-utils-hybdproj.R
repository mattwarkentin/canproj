test_that("hybrid cuttrend is generated with nagg 1", {
  expect_equal(
    get_hybd_cuttrend(0.1, 10, 1, 0.5),
    c(rep(0.1, 5), 0.6, rep(1.0, 4))
  )

  expect_equal(
    get_hybd_cuttrend(0.1, 5, 1, 0.5),
    c(rep(0.1, 5))
  )
})

test_that("hybrid cuttrend is generated with nagg 2", {
  expect_equal(
    get_hybd_cuttrend(0.2, 10, 2, 0.1),
    c(rep(0.2, 3), 0.4, 0.6, 0.8, rep(1.0, 4))
  )

  expect_equal(
    get_hybd_cuttrend(0.2, 3, 2, 0.5),
    c(rep(0.2, 3))
  )
})

test_that("hybrid cuttrend is generated with nagg 3", {
  expect_equal(
    get_hybd_cuttrend(0.1, 10, 3, 0.2),
    c(rep(0.1, 2), 0.7, rep(1.0, 7))
  )

  expect_equal(
    get_hybd_cuttrend(0.1, 2, 3, 0.5),
    rep(0.1, 2)
  )
})

test_that("hybrid cuttrend is generated with nagg 4", {
  expect_equal(
    get_hybd_cuttrend(0.1, 8, 4, 0.2),
    c(0.1, 0.9, rep(1.0, 6))
  )

  expect_equal(
    get_hybd_cuttrend(0.1, 1, 4, 0.5),
    0.1
  )
})

test_that("hybrid cuttrend is generated with nagg 5", {
  expect_equal(
    get_hybd_cuttrend(0.2, 8, 5, 0.1),
    c(0.2, 0.7, rep(1.0, 6))
  )

  expect_equal(
    get_hybd_cuttrend(0.3, 1, 5, 0.5),
    0.3
  )
})
