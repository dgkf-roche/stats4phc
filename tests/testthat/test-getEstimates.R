outcome <- as.numeric(rep(c(TRUE, FALSE), 20))
set.seed(3423)
rscore <- rnorm(40, 0.5, 0.5)


test_that("checkBinsQuantiles works", {
  expect_error(
    checkBinsQuantiles(bins = 5, quantiles = 5, score = rscore),
    "together"
  )
  expect_equal(
    checkBinsQuantiles(bins = c(-1.4, 0, 0.5, 1.7), quantiles = NULL, score = rscore),
    list(bins = c(-1.4, 0, 0.5, 1.7), quantiles = NULL)
  )
  expect_error(
    checkBinsQuantiles(bins = c(0.2, 0.5, 1), quantiles = NULL, score = rscore),
    "must be <= min"
  )
  expect_error(
    checkBinsQuantiles(bins = c(-1.4, 0.5, 0.8), quantiles = NULL, score = rscore),
    "must be >= max"
  )
  expect_error(
    checkBinsQuantiles(bins = 0.5, quantiles = NULL, score = rscore),
    "> 1"
  )
  expect_warning(
    checkBinsQuantiles(bins = 5, quantiles = NULL, score = rep(c(1, 2), each = 5)),
    "> the number of unique score values"
  )
  expect_warning(
    checkBinsQuantiles(
      bins = c(0.5, 1.2, 1.5, 2.5),
      quantiles = NULL,
      score = rep(c(1, 2), each = 5)
    ),
    "> the number of unique score values"
  )
  expect_warning(
    checkBinsQuantiles(bins = NULL, quantiles = 5, score = rep(c(1, 2), each = 5)),
    "Using the quantile method"
  )
})


test_that("getBINNEDest adds outcome.low, outcome.high, and midquantile", {
  res <- getBINNEDest(outcome, rscore, quantiles = 10)
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 10)

  # test quantiles
  expect_snapshot(res)

  # test bins
  res <- getBINNEDest(outcome, rscore, bins = 8)
  expect_snapshot(res)

  expect_null(attr(res, "errorbar"))
})



test_that("getERRORest returns errorbar data of correct dimensions", {
  binn <- getBINNEDest(outcome, rscore, quantiles = 10, errorbar.sem = 1)
  err <- attr(binn, "errorbar")
  expect_equal(nrow(err), 10)
  expect_equal(ncol(err), 5)
  expect_snapshot(as.data.frame(err))
})


test_that("getPAVAest returns estimates of correct dimensions", {
  res <- getPAVAest(outcome, rscore)
  expect_equal(nrow(res), length(rscore))
  expect_equal(ncol(res), 4)
  expect_equal(is.data.frame(res), TRUE)
  expect_snapshot(res)

  res <- getPAVAest(outcome, rscore, low_nonevents = 5)
  expect_snapshot(res)
  res <- getPAVAest(outcome, rscore, high_events = 5)
  expect_snapshot(res)
  res <- getPAVAest(outcome, rscore, hilo_obs = 5)
  expect_snapshot(res)
})

test_that("getGAMest returns estimates of correct dimensions", {
  
  res1 <- getGAMest(outcome, rscore)
  expect_equal(nrow(res1), length(rscore))
  expect_equal(ncol(res1), 4)

  res2 <- getGAMest(outcome, rscore, fitonPerc = FALSE)
  res3 <- getGAMest(outcome, rscore, logscores = TRUE)
  res4 <- getGAMest(outcome[rscore > 0], rscore[rscore > 0], logscores = TRUE, fitonPerc = FALSE)

  expect_snapshot(res1)
  expect_snapshot(res2)
  expect_snapshot(res3)
  expect_snapshot(res4)

})

test_that("getCGAMest works", {
  
  res1 <- getCGAMest(outcome, rscore)
  expect_equal(nrow(res1), length(rscore))
  expect_equal(ncol(res1), 4)

  res2 <- getCGAMest(outcome, rscore, fitonPerc = FALSE)
  
  expect_snapshot(res1)
  expect_snapshot(res2)
})


test_that("getMSPLINEest returns estimates of correct dimensions", {
  
  res1 <- getMSPLINEest(outcome, rscore)
  
  expect_equal(nrow(res1), length(rscore))
  expect_equal(ncol(res1), 4)

  res2 <- getMSPLINEest(outcome, rscore, fitonPerc = FALSE)
  
  expect_snapshot(res1)
  expect_snapshot(res2)
})


test_that("getEsts returns estimates of correct dimensions and correct number of methods", {
  
  expect_error(
    getEsts(
      methods = list(my_est = function(outcome, score) {"dummy"}),
      outcome = outcome,
      score = rscore
    ),
    "data.frame of 4 columns"
  )
  
  expect_error(
    getEsts(
      methods = list(my_est = function(outcome, score) {
        data.frame(score = "a", percentile = "b", outcome = "c", estimate = 0)
      }),
      outcome = outcome,
      score = rscore
    ),
    "not numeric"
  )
  
  expect_error(
    getEsts(methods = methodCheck("dummy"), outcome = outcome, score = score),
    "is not yet available"
  )
  
  expect_error(
    getEsts(
      methods = methodCheck(list(dummy = list(method = "dummy"))),
      outcome = outcome, score = rscore
    ),
    "is not yet available"
  )
  
  expect_error(
    getEsts(
      methods = methodCheck(list(user_gam = function(x) {})),
      outcome = outcome, score = rscore
    ),
    "exactly two arguments"
  )
  
  expect_error(
    getEsts(
      methods = methodCheck(list(asis1 = list(method = "asis"), asis2 = list(method = "asis"))),
      outcome = outcome, score = rscore
    ),
    "just once"
  )
  
  res <- getEsts(methods = list(GAM = list(method = "gam")), outcome = outcome, score = rscore)

  expect_equal(nrow(res$plotdata), length(rscore))
  expect_equal(ncol(res$plotdata), 5)
  expect_equal(length(names(table(res$plotdata$method))), 1)
  expect_snapshot(res$plotdata)

  res <- getEsts(
    methods = list("GAM" = list(method = "gam"), "PAVA" = list(method = "pava")),
    outcome = outcome,
    score = rscore
  )
  expect_equal(length(names(table(res$plotdata$method))), 2)
  expect_snapshot(res$plotdata)
  expect_equal(res$idx.step, c("GAM" = F, "PAVA" = T))
  expect_equal(res$idx.asis, c("GAM" = F, "PAVA" = F))

  res <- getEsts(
    methods = list(
      "GAM" = list(method = "gam"),
      "PAVA" = list(method = "pava"),
      "MSPLINE" = list(method = "mspline")
    ),
    outcome = outcome,
    score = rscore
  )
  expect_equal(length(names(table(res$plotdata$method))), 3)
  expect_snapshot(res$plotdata)

  res <- getEsts(
    methods = list(
      "GAM" = list(method = "gam"),
      "PAVA" = list(method = "pava"),
      "MSPLINE" = list(method = "mspline"),
      "BINNED" = list(method = "binned")
    ),
    outcome = outcome,
    score = rscore
  )
  expect_equal(length(names(table(res$plotdata$method))), 4)
  expect_snapshot(res$plotdata)
  expect_null(res$errorbardata)
  expect_equal(res$idx.step, c("GAM" = F, "PAVA" = T, "MSPLINE" = F, "BINNED" = T))
  expect_equal(res$idx.asis, c("GAM" = F, "PAVA" = F, "MSPLINE" = F, "BINNED" = F))
  expect_equal(res$idx.binned, c("GAM" = F, "PAVA" = F, "MSPLINE" = F, "BINNED" = T))
  expect_equal(res$idx.pava, c("GAM" = F, "PAVA" = T, "MSPLINE" = F, "BINNED" = F))

  res <- getEsts(
    methods = list(
      "GAM" = list(method = "gam"),
      "BINNED" = list(method = "binned", errorbar.sem = 0.87)
    ),
    outcome = outcome,
    score = rscore
  )
  expect_snapshot(as.data.frame(res$errorbardata))

  res <- getEsts(
    methods = list(
      "g" = list(method = "gam"),
      "b" = list(method = "binned"),
      "r" = list(method = "asis")
    ),
    outcome = outcome,
    score = rscore
  )
  expect_snapshot(as.data.frame(res$plotdata))
  expect_equal(res$idx.step, c("g" = F, "b" = T, "r" = F))
  expect_equal(res$idx.asis, c("g" = F, "b" = F, "r" = T))
  expect_equal(res$idx.binned, c("g" = F, "b" = T, "r" = F))
  expect_equal(res$idx.pava, c("g" = F, "b" = F, "r" = F))

})

test_that("asis method works in getEsts", {
  res <- getEsts(
    methods = methodCheck("asis"),
    outcome = outcome,
    score = rscore
  )
  expect_snapshot(as.data.frame(res$plotdata))
})


test_that("getConstraints works", {
  out <- getConstraints(outcome = outcome, rscore = rscore, low_events = 3, high_nonevents = 3)
  expect_snapshot(out)

  out <- getConstraints(
    outcome = outcome, rscore = rscore,
    low_events = 0, low_nonevents = 0,
    high_events = 0, high_nonevents = 0,
    hilo_obs = 0
  )
  expect_snapshot(out)

  expect_warning(
    getConstraints(
      outcome = outcome, rscore = rscore,
      low_events = 1, low_nonevents = 1
    )
  )
  out <- suppressWarnings(
    getConstraints(
      outcome = outcome, rscore = rscore,
      low_events = 1, low_nonevents = 1
    )
  )
  expect_snapshot(out)


  expect_warning(
    getConstraints(
      outcome = outcome, rscore = rscore,
      high_events = 1, high_nonevents = 1
    )
  )
  out <- suppressWarnings(
    getConstraints(
      outcome = outcome, rscore = rscore,
      high_events = 1, high_nonevents = 1
    )
  )
  expect_snapshot(out)
})



test_that("summaryTable argument returns appropriately sized dataframe when requested.", {
  
  expect_error(
    getSummaries(outcome, rscore, bins = numeric(0), quantiles = NULL),
    "Unrecognized option"
  )
  
  expect_equal(nrow(getSummaries(outcome, rscore, quantiles = 10)$binlvl), 10)
  expect_equal(ncol(getSummaries(outcome, rscore, quantiles = 10)$binlvl), 9)

  expect_equal(nrow(getSummaries(outcome, rscore, bins = 10)$binlvl), 10)
  expect_equal(ncol(getSummaries(outcome, rscore, bins = 10)$binlvl), 9)
  expect_error(getSummaries(outcome, rscore, quantiles = 0, bins = 0))

  expect_equal(is.data.frame(getSummaries(outcome, rscore, quantiles = 10)$obslvl), TRUE)
  expect_equal(is.data.frame(getSummaries(outcome, rscore, quantiles = 10)$binlvl), TRUE)

  # test quantiles
  out <- getSummaries(outcome, rscore, quantiles = 7)
  expect_snapshot(out$obslvl)
  expect_snapshot(as.data.frame(out$binlvl))

  # test bins
  out <- getSummaries(outcome, rscore, bins = 8)
  expect_snapshot(out$obslvl)
  expect_snapshot(as.data.frame(out$binlvl))

  # test bins - 2
  out <- getSummaries(outcome, rscore, bins = c(0, 0.4, 0.75, 1))
  expect_snapshot(out$obslvl)
  expect_snapshot(as.data.frame(out$binlvl))

  # test right
  out <- getSummaries(outcome, rscore, bins = 8, right = FALSE)
  expect_snapshot(out$obslvl)
  expect_snapshot(as.data.frame(out$binlvl))

})
