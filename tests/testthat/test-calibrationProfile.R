outcome <- rep(c(TRUE, FALSE), 20)
set.seed(3423)
rscore <- rnorm(40, 0.5, 0.5)

test_that("calibrationProfile works", {
  
  expect_error(
    calibrationProfile(outcome = outcome, score = rscore, methods = "asis"),
    "not suitable"
  )
  
  out1 <- calibrationProfile(outcome = outcome, score = rscore)
  out2 <- calibrationProfile(outcome = outcome, score = rscore, plot.raw = FALSE)
  out3 <- calibrationProfile(outcome = outcome, score = rscore, rev.order = TRUE)
  
  out4 <- calibrationProfile(
    outcome = outcome, score = rscore,
    margin.type = "histogram", bins = 50
  )
  
  out5 <- calibrationProfile(
    outcome = outcome, score = rscore,
    include = c("citl", "datapoints", "rug", "loess")
  )
  
  out6 <- calibrationProfile(
    outcome = outcome, score = rscore,
    methods = c("BINNED", "PAVA", "MSPLINE", "gam")
  )
  
  out7 <- calibrationProfile(
    outcome = outcome, score = rscore,
    methods = "gam", include = NULL
  )
  
  expect_snapshot(as.data.frame(out1$data))
  expect_snapshot(out1$citl)

  expect_snapshot(as.data.frame(out2$data))
  expect_snapshot(out2$citl)

  expect_snapshot(as.data.frame(out3$data))
  expect_snapshot(out3$citl)

  expect_s3_class(out4$plot, "gtable")

  expect_snapshot(as.data.frame(out5$data))
  expect_snapshot(out5$citl)

  expect_snapshot(as.data.frame(out6$data))
  expect_snapshot(out6$citl)

  expect_snapshot(as.data.frame(out7$data))
  expect_snapshot(out7$citl)

  # expect_snapshot_file(save_png(out1$plot), "p1.png")
  # expect_snapshot_file(save_png(out2$plot), "p2.png")
  # expect_snapshot_file(save_png(out3$plot), "p3.png")
  # expect_snapshot_file(save_png(out5$plot), "p5.png")
  # expect_snapshot_file(save_png(out6$plot), "p6.png")
  # expect_snapshot_file(save_png(out7$plot), "p7.png")
})

test_that("calibrationProfile with user defined function works", {
  p1 <- calibrationProfile(outcome = outcome, score = rscore, methods = "GAM")
  p2 <- calibrationProfile(
    outcome = outcome, score = rscore,
    methods = list(user_gam = function(outcome, score) {
      perc <- ecdf(score)(score)
      s <- mgcv::s
      m <- mgcv::gam(
        outcome ~ s(perc, bs = "tp"),
        family = "binomial",
        method = "REML"
      )
      preds <- mgcv::predict.gam(m, type = "response")
      data.frame(score = score, percentile = perc, outcome = outcome, estimate = preds)
    })
  )
  p2$data$method <- "GAM"

  expect_equal(p1$data, p2$data)

  # expect_snapshot_file(save_png(p2$plot), "udf.png")
})
