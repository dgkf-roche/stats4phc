outcome <- rep(c(TRUE, FALSE), 20)
set.seed(3423)
rscore <- rnorm(40, 0.5, 0.5)

test_that("riskProfile works", {
  
  out1 <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = c("MSPLINE", "GAM", "BINNED", "PAVA"),
    rev.order = FALSE
  )
  
  out2 <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = c("MSPLINE", "GAM", "BINNED", "PAVA"),
    rev.order = TRUE
  )
  
  out3 <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = c("GAM", "ASIS"),
    rev.order = FALSE,
    prev.adj = 0.37
  )
  
  out4 <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = "asis",
    show.best.pv = TRUE,
    show.nonparam.pv = FALSE,
    include = c("PPV", "NPV", "1-NPV")
  )
  
  expect_snapshot(as.data.frame(out1$data))
  expect_null(out1$errorbar)

  expect_snapshot(as.data.frame(out2$data))
  expect_null(out2$errorbar)

  expect_snapshot(as.data.frame(out3$data))
  expect_null(out3$errorbar)

  expect_equal(nrow(out4$data), 246)
  expect_null(out4$errorbar)

  # expect_snapshot_file(save_png(out1$plot), "riskProfile.png")
  # expect_snapshot_file(save_png(out2$plot), "riskProfile_rev_order.png")
  # expect_snapshot_file(save_png(out3$plot), "riskProfile_prev_best.png")
  # expect_snapshot_file(save_png(out4$plot), "riskProfile-empty.png")
})


test_that("BINNED method with errorbar.sem and passing method parameters works", {
  out <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = list(bin = list(method = "BINNED", errorbar.sem = 1.2))
  )

  expect_snapshot(as.data.frame(out$data))
  expect_snapshot(as.data.frame(out$errorbar))

  # expect_snapshot_file(save_png(out$plot), "errorbar.png")
})

test_that("plot.raw and rev.order works", {
  p1 <- riskProfile(outcome = outcome, score = rscore, plot.raw = F, rev.order = F)
  p2 <- riskProfile(outcome = outcome, score = rscore, plot.raw = T, rev.order = F)
  p3 <- riskProfile(outcome = outcome, score = rscore, plot.raw = F, rev.order = T)
  p4 <- riskProfile(outcome = outcome, score = rscore, plot.raw = T, rev.order = T)

  expect_snapshot(as.data.frame(p1$data))

  expect_snapshot(as.data.frame(p2$data))

  expect_snapshot(as.data.frame(p3$data))

  expect_snapshot(as.data.frame(p4$data))

  # expect_snapshot_file(save_png(p1$plot), "p1.png")
  # expect_snapshot_file(save_png(p2$plot), "p2.png")
  # expect_snapshot_file(save_png(p3$plot), "p3.png")
  # expect_snapshot_file(save_png(p4$plot), "p4.png")
})


test_that("User defined estimation functions work", {
  p1 <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = "GAM",
    show.nonparam.pv = FALSE,
    show.best.pv = FALSE
  )

  p2 <- riskProfile(
    outcome = outcome,
    score = rscore,
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
    }),
    show.nonparam.pv = FALSE,
    show.best.pv = FALSE
  )
  p2$data$method <- "GAM"

  expect_equal(p1$data, p2$data)

  # expect_snapshot_file(save_png(p2$plot), "udf.png")
})

test_that("CGAM works", {
  out <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = "CGAM"
  )
  expect_snapshot(as.data.frame(out$data))

  # expect_snapshot_file(save_png(out$plot), "cgam.png")
})

test_that("comparing single predictive quantity works", {
  out <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = list(
      "gam_3" = list(method = "gam", k = 3),
      "gam_4" = list(method = "gam", k = 4),
      "gam_7" = list(method = "gam", k = 7)
    ),
    show.best.pv = FALSE,
    show.nonparam.pv = FALSE,
    include = "PC"
  )
  expect_snapshot(as.data.frame(out$data))

  # expect_snapshot_file(save_png(out$plot), "compare_pc.png")
})

test_that("single PV (not PC) works", {
  out <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = "gam",
    show.best.pv = TRUE,
    include = "PPV"
  )
  expect_snapshot(as.data.frame(out$data))

  # expect_snapshot_file(save_png(out$plot), "single_pv.png")
})

test_that("arguments show.best.pv and show.nonparam.pv works", {
  out <- riskProfile(
    outcome = outcome,
    score = rscore,
    methods = "gam",
    show.best.pv = FALSE,
    show.nonparam.pv = FALSE,
    include = c("PPV", "NPV")
  )
  expect_snapshot(as.data.frame(out$data))

  # expect_snapshot_file(save_png(out$plot), "show_arguments.png")
})
