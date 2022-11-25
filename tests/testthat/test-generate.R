test_that("The dimensions are numeric", {
  expect_error(generate(L = "a", K = 1))
  expect_error(generate(L = 1, K = "a"))
})

test_that("The parameters are well defined", {
  expect_error(generate(L = 1, K = 1, g_endo = -1))
  expect_error(generate(L = 1, K = 1, g_endo = 2))
})
