test_that("zipFastener checks inputs", {
  expect_error(
    zipFastener(data.frame(x = 1:4), data.frame(x = 1:3)),
    "the no. of rows"
  )
  expect_error(
    zipFastener(data.frame(x = 1:4, y = 1:4), data.frame(x = 1:4)),
    "the no. of columns"
  )
})

# other tests done within riskProfile function
