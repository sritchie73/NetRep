context("Permutation Testing")

n <- 100000
normData <- rnorm(n=n)

test_that("perm.test approximates normal for n=100000", {
  expect_equal(
    perm.test(normData, -1.644854, 1, n, F, "less"), 
    pnorm(-1.644854, lower.tail=TRUE),
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, n, F, "greater"), 
    pnorm(1.644854, lower.tail=FALSE),
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, n, F, "two.sided"), 
    pnorm(1.644854, lower.tail=FALSE)*2,
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, n, F, "two.sided"), 
    pnorm(-1.644854, lower.tail=TRUE)*2,
    tolerance = 0.05
  )
})

