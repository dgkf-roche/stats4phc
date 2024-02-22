outcome <- rep(c(TRUE, FALSE), 20)
set.seed(5112)
rscore <- rnorm(40, 0.5, 0.5)


test_that("inputCheck works", {
  
  expect_error(inputCheck(c(1, 0, 1), c("a", "b", "c")), "numeric")
  expect_error(inputCheck(c("a", "b"), c(0.1, 0.2)), "numeric")
  expect_error(inputCheck(c(0, 1, 0), c(0.1, 0.2)), "same lengths")
  expect_error(inputCheck(c(1, 0), c(0.1, 0.2)), "observations")

  res <- suppressWarnings(inputCheck(c(1, 0, 1, NA, 0), c(0.1, NA, 0.2, 0.3, 0.4)))
  expect_equal(res, list(outcome = c(1, 1, 0), score = c(0.1, 0.2, 0.4)))

  res <- suppressWarnings(inputCheck(c(1, 0, 1, NA, 0), c(0.1, NA, 0.2, 0.3, 0.4)))
  expect_equal(res, list(outcome = c(1, 1, 0), score = c(0.1, 0.2, 0.4)))
  
  expect_error(inputCheck(c(1, 0, 2), c(0.1, 0.2, 0.3)), "as binary")
  
  expect_warning(
    inputCheck(c(1, rep(0, 49)), 1:50),
    "There is a low frequency"
  )
  
  res <- inputCheck(outcome, rscore)

  expect_equal(length(res$outcome), length(res$score))
  expect_snapshot(res)
  expect_equal(is.numeric(res$outcome), TRUE)
  expect_equal(is.numeric(res$score), TRUE)

  expect_equal(
    inputCheck(rep(c(TRUE, FALSE), 5), 1:10 / 10),
    list(outcome = as.numeric(rep(c(TRUE, FALSE), 5)), score = 1:10 / 10)
  )

  rscore[1] <- NA
  expect_warning(inputCheck(outcome, rscore))
})



test_that("methods returns appropriate length list of lists", {
  expect_equal(is.list(methodCheck(methods = c("GAM"))), TRUE)
  expect_equal(length(methodCheck(methods = c("GAM", "PAVA", "asis"))), 3)
  expect_equal(is.list(methodCheck(methods = list("GAM" = list(method = "GAM")))), TRUE)
  expect_equal(
    length(
      methodCheck(methods = list("GAM" = list(method = "GAM"), "PAVA" = list(method = "PAVA")))
    ),
    2
  )

  ms <- c("asis", "GAM", "MSPLINE", "BINNED", "PAVA")
  expect_equal(
    methodCheck(ms),
    structure(as.list(tolower(ms)), names = ms)
  )

  expect_equal(
    methodCheck(list(gam = list(method = "gam", k = 5))),
    list(gam = list(method = "gam", k = 5))
  )

  expect_error(methodCheck(NULL))
  expect_equal(methodCheck(c("GAM", "")), list(GAM = "gam"))
  expect_error(methodCheck(1))

  expect_error(methodCheck(list(GAM = c())), "as a list or a function")

  expect_error(methodCheck(c("GAM", "GAM")), "unique")
  expect_error(methodCheck(list(GAM = list(), GAM = list())), "uniquely")

  expect_error(methodCheck(list(m1 = list("gam"))), "named")
  expect_error(
    methodCheck(list(m1 = list(method = 1))),
    "specifying one of the predefined estimation functions"
  )
  expect_error(
    methodCheck(list(m1 = list(method = NULL))),
    "specifying one of the predefined estimation functions"
  )
  
  expect_error(methodCheck(c("gamm")), "not yet available")
})



test_that("orderInputs works as expected", {
  expect_snapshot(orderInputs(outcome = outcome, score = rscore, rev.order = FALSE))
  expect_snapshot(orderInputs(outcome = outcome, score = rscore, rev.order = TRUE))
})


test_that("add0thPercPC works", {
  df <- data.frame(
    percentile = c(0.1, 0.5, 1, 0.2, 0.4, 0.6),
    estimate = c(6, 5, 4, 5, 5, 5),
    method = rep(c("m1", "m2"), each = 3),
    pv = "pv"
  )
  out <- add0thPercPC(df)
  expect_equal(
    dplyr::bind_rows(
      data.frame(percentile = 0, estimate = 6, method = "m1", pv = "PC", score = NA, outcome = NA),
      df[1:3, ],
      data.frame(percentile = 0, estimate = 5, method = "m2", pv = "PC", score = NA, outcome = NA),
      df[4:6, ]
    ),
    out
  )
})


test_that("add0thPercPV works", {
  df1 <- data.frame(
    percentile = c(0.1, 0.5, 1, 0.2, 0.4, 0.6),
    pvValue = c(6, 5, 4, 5, 5, 5),
    method = rep(c("m1", "m2"), each = 3),
    pv = "pv"
  )
  df2 <- df1
  df2$pv <- "pv2"
  df <- dplyr::bind_rows(df1, df2)
  out <- add0thPercPV(df)

  comp1 <- dplyr::bind_rows(
    data.frame(percentile = 0, pvValue = 6, method = "m1", pv = "pv", score = NA, estimate = NA),
    df1[1:3, ],
    data.frame(percentile = 0, pvValue = 6, method = "m1", pv = "pv2", score = NA, estimate = NA),
    df2[1:3, ]
  )

  comp2 <- dplyr::bind_rows(
    data.frame(percentile = 0, pvValue = 5, method = "m2", pv = "pv", score = NA, estimate = NA),
    df1[4:6, ],
    data.frame(percentile = 0, pvValue = 5, method = "m2", pv = "pv2", score = NA, estimate = NA),
    df2[4:6, ]
  )

  expect_equal(bind_rows(comp1, comp2), out)
})


test_that("listMethodCheck works", {
  
  expect_error(listMethodCheck(list(method = "my_method")), "not yet available")
  
  expect_error(listMethodCheck(function(x) {x}), "exactly two arguments")
  
})
