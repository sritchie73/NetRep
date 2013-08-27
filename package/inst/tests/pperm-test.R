context("Functions on Distributions from Permutation Testing")

normData <- rnorm(n=10000)
small <- sort(runif(4))
approx <- function(...) { equals(..., tolerance=0.1)}

test_that("pperm approximates normal for n=10000", {
  expect_that(pperm(normData, -1.644854), approx(pnorm(-1.644854)))
  expect_that(pperm(normData, 1.644854), approx(pnorm(1.644854)))
  expect_that(pperm(normData, -1.644854, lower.tail=FALSE), approx(pnorm(1.644854)))  
  expect_that(pperm(normData, 1.644854, lower.tail=FALSE), approx(pnorm(-1.644854)))  
  expect_that(pperm(normData, -1.644854, log.p=TRUE), approx(pnorm(-1.644854, log.p=TRUE)))
})

test_that("qperm returns similar pvalues to pnorm for n=10000", {
  expect_that(qperm(normData, 0.05), approx(qnorm(0.05)))
  expect_that(qperm(normData, 0.01), approx(qnorm(0.01)))
  expect_that(qperm(normData, 0.95), approx(qnorm(0.95)))
  expect_that(qperm(normData, -0.05129329, log.p=TRUE), approx(qnorm(-0.05129329, log.p=TRUE)))
  expect_that(qperm(normData, -2.995732, log.p=TRUE), approx(qnorm(-2.995732, log.p=TRUE)))
})

test_that("qperm and pperm agree on p-values calculation", {
  ps <- sapply(small, function(x) { pperm(small, x) })
  qs <- sapply(ps, function(x) { qperm(small, x) })

  expect_identical(small, qs)
})

test_that("qperm for range errors", {
  expect_identical(qperm(small, 0.1), small[1])
  expect_identical(qperm(small, 0.9), small[4])
  expect_identical(qperm(small, log(0.1), log.p=TRUE), small[1])
  expect_identical(qperm(small, log(0.9), log.p=TRUE), small[4])
})

test_that("qperm handles multiple returns correctly", {  
  expect_warning(qperm(small, 0.5))
})
