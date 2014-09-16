context("Functions on Distributions from Permutation Testing")

normData <- rnorm(n=10000)
small <- sort(runif(4))
# testing if approximately the same value from a normal
normApprox <- function(...) { equals(..., tolerance = 0.1) }
# testing if approximately the same p-value
pApprox <- function(...) { equals(..., tolerance = 0.05) }

test_that("pperm approximates normal for n=10000", {
  expect_that(pperm(normData, -1.644854, 1), pApprox(pnorm(-1.644854, lower.tail=FALSE)))
  expect_that(pperm(normData, 1.644854, 1), pApprox(pnorm(1.644854, lower.tail=FALSE)))  
})

