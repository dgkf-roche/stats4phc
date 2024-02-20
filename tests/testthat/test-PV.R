outcome <- as.numeric(rep(c(TRUE, FALSE), 20))
set.seed(3423)
rscore <- rnorm(40, 0.5, 0.5)
cdf.fit <- ecdf(rscore)
percentile <- cdf.fit(rscore)


test_that("getPVdata works", {
  out <- getPVdata(
    outcome = outcome,
    score = rscore,
    methods = list(BINNED = list(method = "binned"), MSPLINE = list(method = "mspline"))
  )
  expect_snapshot(out)

  out <- getPVdata(
    outcome = outcome,
    score = rscore,
    methods = list(
      BINNED = list(method = "binned"), MSPLINE = list(method = "mspline"),
      pv = list(method = "pava"), bin2 = list(method = "binned", bins = 8)
    )
  )
  expect_snapshot(out)
})

test_that("getPVdata works with asis", {
  out <- getPVdata(
    outcome = as.numeric(outcome),
    score = rscore,
    methods = methodCheck("asis")
  )
  expect_snapshot(out)
})


test_that("nonParametricPV works", {
  pvdata <- nonParametricPV(score = rscore, as.numeric(outcome))

  expect_equal(nrow(pvdata), length(rscore))
  expect_equal(ncol(pvdata), 5)
  expect_equal(all(pvdata$score == rscore), TRUE)
  expect_equal(all(pvdata$Percentile == percentile), TRUE)

  expect_snapshot(pvdata)
})


test_that("adjust prevalence works", {
  cdf1 <- ecdf(rscore[outcome == 1])
  cdf2 <- ecdf(rscore[outcome == 0])

  ap <- adjPrevPerc(perc = percentile, prev.new = 0.8, cdf.control = cdf1, cdf.case = cdf2)
  expect_equal(head(ap, 6), c(0.8, 0.3, 0.47, 0.4, 0.62, 0.26))
  expect_equal(tail(ap, 6), c(0.78, 0.78, 0.62, 0.3, 0.61, 0.41))

  ae <- round(adjPrevEst(x = 1:10, prev = 0.2, prev.new = 0.7), 2)
  expect_equal(ae, c(1, 1.06, 1.08, 1.09, 1.09, 1.1, 1.1, 1.1, 1.11, 1.11))

  ae <- round(adjPrevEst(x = 1:10, prev = 0.8, prev.new = 0.3), 2)
  expect_equal(ae, c(1, -0.27, -0.19, -0.17, -0.15, -0.15, -0.14, -0.14, -0.14, -0.14))
})


test_that("best ppv/mnpv works", {
  ppv <- round(bestPPV(perc = percentile, prev = mean(outcome)), 2)
  expect_equal(head(ppv, 6), c(1, 0.57, 0.77, 0.67, 1, 0.51))
  expect_equal(tail(ppv, 6), c(1, 1, 1, 0.56, 1, 0.69))

  mnpv <- round(bestMNPV(perc = percentile, prev = mean(outcome)), 2)
  expect_equal(head(mnpv, 6), c(0.5, 0, 0, 0, 0.23, 0))
  expect_equal(tail(mnpv, 6), c(0.41, 0.43, 0.29, 0, 0.2, 0))
})


test_that("best sens/spec works", {
  sn <- bestSens(perc = percentile, prev = mean(outcome))
  expect_equal(head(sn, 6), c(0, 1, 1, 1, 0.7, 1))
  expect_equal(tail(sn, 6), c(0.3, 0.25, 0.6, 1, 0.75, 1))

  sp <- bestSpec(perc = percentile, prev = mean(outcome))
  expect_equal(head(sp, 6), c(1, 0.25, 0.7, 0.5, 1, 0.05))
  expect_equal(tail(sp, 6), c(1, 1, 1, 0.2, 1, 0.55))
})
