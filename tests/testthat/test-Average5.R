test_that("ave5proj works when aggregating data", {
  cdat <- data.frame(matrix(1:190, nrow = 19, ncol = 10))
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 15))
  startp <- 2000

  projection <- data.frame(matrix(
    c(
      seq(from = 10, to = 1900, by = 10),
      rep(seq(from = 1245, to = 1425, by = 10), 5)
    ),
    nrow = 19,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  expected <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = NULL,
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
      rep(seq(from = 1150, to = 1330, by = 10), 5)
    ),
    nrow = 19,
    ncol = 15
  ))
  colnames(projection) <- c(1990:2004)
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

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
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

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
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

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

test_that("ave5prog gets projections with default standard population", {
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
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  ave5_object <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )

  expect_equal(ave5proj.getproj(pdat, ave5_object), projection)
})

test_that("ave5proj gets projections with specified standard population", {
  pdat <- data.frame(matrix(10000, nrow = 19, ncol = 5))

  projection <- data.frame(matrix(
    seq(from = 10, to = 480, by = 5),
    nrow = 19,
    ncol = 5
  ))
  colnames(projection) <- c(1990:1994)
  rownames(projection) <- c(
    paste0(seq(0, 85, by = 5), "-", seq(4, 89, by = 5)),
    "90+"
  )

  ave5_object <- list(
    agsproj = projection,
    startp = 2000,
    sum5 = TRUE,
    noypred = 5
  )

  stdpop <- c(rep(0.08, 9), rep(0.025, 10))

  out <- matrix(
    c(
      c(40.975, 133.125, 225.275, 317.425, 409.575),
      c(104, 285, 466, 646, 826)
    ),
    nrow = 5,
    ncol = 2
  )
  rownames(out) <- 1990:1994
  colnames(out) <- c("asr", "case")

  expect_equal(ave5proj.getproj(pdat, ave5_object, stdpop), out)
})

test_that("ave5proj summary accepts only ave5proj class objects", {
  m <- matrix(1, nrow = 15, ncol = 5)

  expect_error(
    summary.ave5proj(m),
    "Variable \"ave5proj.object\" must be of type \"ave5proj\""
  )
})
