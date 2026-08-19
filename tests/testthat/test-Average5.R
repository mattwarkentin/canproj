test_that("ave5proj works when aggregating data", {
  cdat <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 15))
  startp <- 2000

  projection <- data.frame(matrix(
    c(
      seq(from = 10, to = 1900, by = 10),
      rep(seq(from = 1340, to = 1520, by = 10), 5)
    ),
    nrow = 19,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- 1:19

  expected <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )
  class(expected) <- "ave5proj"
  attr(expected, "Call") <- as.call(str2lang("ave5proj(cdat, pdat, startp)"))

  expect_equal(ave5proj(cdat, pdat, startp), expected)
})

test_that("ave5proj works when averaging data", {
  cdat <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 15))
  startp <- 2000

  projection <- data.frame(matrix(
    c(
      seq(from = 10, to = 1900, by = 10),
      rep(seq(from = 1340, to = 1520, by = 10), 5)
    ),
    nrow = 19,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- 1:19

  expected <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )
  class(expected) <- "ave5proj"
  attr(expected, "Call") <- as.call(str2lang(
    "ave5proj(cdat, pdat, startp, sum5 = T)"
  ))

  expect_equal(ave5proj(cdat, pdat, startp, sum5 = T), expected)
})

test_that("ave5proj works with non-standard age groups", {
  cdat <- data.frame(matrix(1:100, nrow = 10, ncol = 10))
  pdat <- data.frame(matrix(10000, nrow = 10, ncol = 15))
  startp <- 2000

  out <- ave5proj(cdat, pdat, startp)

  projection <- data.frame(matrix(
    c(
      seq(from = 10, to = 1000, by = 10),
      rep(seq(from = 710, to = 800, by = 10), 5)
    ),
    nrow = 10,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- 1:10

  expect_equal(out$agsproj, projection)
})

test_that("ave5proj creates summary with printcall", {
  cdat <- data.frame(matrix(1:38, nrow = 19, ncol = 2))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 4))
  startp <- 2000

  projection <- data.frame(matrix(
    seq(from = 10, to = 760, by = 10),
    nrow = 19,
    ncol = 4
  ))
  colnames(projection) <- c(1990:1993)
  rownames(projection) <- 1:19

  ave5_obj <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = NULL,
    noypred = 2
  )
  class(ave5_obj) <- "ave5proj"
  attr(ave5_obj, "Call") <- as.call(str2lang("ave5proj(cdat, pdat, startp)"))

  out <- capture.output(summary.ave5proj(ave5_obj, printcall = TRUE))

  expect_match(out, "Prediction done with:", all = FALSE)
  expect_match(out, "Method:\\s+Five-Year Average", all = FALSE)
  expect_match(out, "Age-Specific Rate by:\\s+5-year period", all = FALSE)
  expect_match(out, "Number of prediction years:\\s+2", all = FALSE)
  expect_match(out, "Call:\\s+ave5proj\\(cdat, pdat, startp\\)", all = FALSE)
})

test_that("ave5proj creates summary with averaged yearly rates", {
  cdat <- data.frame(matrix(1:38, nrow = 19, ncol = 2))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 4))
  startp <- 2000

  projection <- data.frame(matrix(
    seq(from = 10, to = 760, by = 10),
    nrow = 19,
    ncol = 4
  ))
  colnames(projection) <- c(1990:1993)
  rownames(projection) <- 1:19

  ave5_obj <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = T,
    noypred = 2
  )
  class(ave5_obj) <- "ave5proj"
  attr(ave5_obj, "Call") <- as.call(str2lang("ave5proj(cdat, pdat, startp)"))

  out <- capture.output(summary.ave5proj(ave5_obj, printcall = TRUE))

  expect_match(out, "Prediction done with:", all = FALSE)
  expect_match(out, "Method:\\s+Five-Year Average", all = FALSE)
  expect_match(
    out,
    "Age-Specific Rate by:\\s+average yearly-rates",
    all = FALSE
  )
  expect_match(out, "Number of prediction years:\\s+2", all = FALSE)
  expect_match(out, "Call:\\s+ave5proj\\(cdat, pdat, startp\\)", all = FALSE)
})

test_that("ave5proj gets projections with default standard population", {
  cdat <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 15))

  projection <- data.frame(matrix(
    c(
      seq(from = 10, to = 1900, by = 10),
      rep(seq(from = 1150, to = 1330, by = 10), 5)
    ),
    nrow = 19,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- 1:19

  ave5_object <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )
  class(ave5_object) <- "ave5proj"

  expect_equal(ave5proj_get_projections(pdat, ave5_object), projection)
})

test_that("ave5proj gets projections with specified standard population", {
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 5))

  projection <- data.frame(matrix(
    seq(from = 10, to = 480, by = 5),
    nrow = 19,
    ncol = 5
  ))
  colnames(projection) <- c(1990:1994)
  rownames(projection) <- 1:19

  ave5_object <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )
  class(ave5_object) <- "ave5proj"

  out <- matrix(
    c(
      c(49.120765, 144.120765, 239.120765, 334.120765, 429.120765),
      c(104, 285, 466, 646, 826)
    ),
    nrow = 5,
    ncol = 2
  )
  rownames(out) <- 1990:1994
  colnames(out) <- c("asr", "case")

  expect_equal(
    ave5proj_get_projections(pdat, ave5_object, stdpop_Canada_2021),
    out
  )
})

test_that("ave5proj summary accepts only ave5proj class objects", {
  m <- matrix(1, nrow = 15, ncol = 5)

  expect_error(
    summary.ave5proj(m),
    "Variable \"ave5proj.object\" must be of type \"ave5proj\""
  )
})

test_that("ave5proj plotting works", {
  cdat <- data.frame((row(matrix(1, 19, 15)) - col(matrix(1, 19, 15)) + 50))
  pdat <- data.frame(matrix(50000, nrow = 19, ncol = 20))

  out <- ave5proj(cdat, pdat, 2020)

  vdiffr::expect_doppelganger("acproj plot", function() {
    plot(out, pdat, stdpop_Canada_2021)
  })
})
