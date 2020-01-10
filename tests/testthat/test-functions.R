set.seed(1)

test_that("dbetasigmoid works", {
  xg <- seq(-10, 10, length.out = 1001)
  dens <- dbetasigmoid(xg, a = 3, b = 10)
  expect_equal(sum(dens * diff(xg)[1]), 1)
})


test_that("rbetasigmoid works", {
  samples <- rbetasigmoid(n = 1e5, a = 5, b = 2, infl = 3, scale = 10)
  expect_equal(mean(samples), mean_betasigmoid(a = 5, b = 2, infl = 3, scale = 10), tolerance = 0.01)
  expect_equal(var(samples), var_betasigmoid(a = 5, b = 2, infl = 3, scale = 10), tolerance = 0.01)
})


test_that("pbetasigmoid works", {
  samples <- rbetasigmoid(n = 1e5, a = 10, b = 17, infl = -1, scale = 0.3)
  pvals <- pbetasigmoid(samples, a = 10, b = 17, infl = -1, scale = 0.3)
  expect_equal(cor(sort(pvals), sort(runif(1e5))), 1, tolerance = 1e-5)
})


test_that("qbetasigmoid works", {
  res <- expect_warning(qbetasigmoid(p = c(0, 0.5, 1, 2), a = 10, b = 17, infl = -1, scale = 0.3))
  expect_equal(res[1], -Inf)
  expect_equal(res[3], Inf)
  expect_equal(res[4], NaN)
})

