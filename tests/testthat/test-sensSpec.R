outcome <- rep(c(TRUE, FALSE), 20)
set.seed(3423)
rscore <- rnorm(40, 0.5, 0.5)

cdf.fit <- ecdf(rscore)
percentile <- cdf.fit(rscore)
df <- data.frame(outcome, rscore, percentile)

set.seed(541)
randf <- df[sample(nrow(df)), ] # randomize scores

test_that("sensSpec can order scores", {
  res <- sensSpec(df$outcome, df$rscore)
  res2 <- sensSpec(randf$outcome, randf$rscore)

  expect_equal(res$data, res2$data)
  expect_snapshot(as.data.frame(res$data))

  # expect_snapshot_file(save_png(res$plot), "sensSpec.png")
})

test_that("sensSpec works", {
  res2 <- sensSpec(outcome = outcome, score = rscore, rev.order = T, plot.raw = F)
  expect_snapshot(as.data.frame(res2$data))

  res3 <- sensSpec(outcome = outcome, score = rscore, rev.order = F, plot.raw = F)
  expect_snapshot(as.data.frame(res3$data))

  res4 <- sensSpec(outcome = outcome, score = rscore, rev.order = F, plot.raw = T)
  expect_snapshot(as.data.frame(res4$data))

  res5 <- sensSpec(outcome = outcome, score = rscore, rev.order = T, plot.raw = T)
  expect_snapshot(as.data.frame(res5$data))

  res6 <- sensSpec(outcome = outcome, score = rscore, methods = c("asis", "gam", "pava"))
  expect_snapshot(as.data.frame(res6$data))

  # expect_snapshot_file(save_png(res2$plot), "sensSpec2.png")
  # expect_snapshot_file(save_png(res3$plot), "sensSpec3.png")
  # expect_snapshot_file(save_png(res4$plot), "sensSpec4.png")
  # expect_snapshot_file(save_png(res5$plot), "sensSpec5.png")
  # expect_snapshot_file(save_png(res6$plot), "sensSpec6.png")
})
