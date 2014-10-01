context("Permutation Testing")

normData <- rnorm(n=100000)

test_that("perm.test approximates normal for n=10000", {
  expect_equal(
    perm.test(normData, -1.644854, 1, "less"), 
    pnorm(-1.644854, lower.tail=TRUE),
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, "greater"), 
    pnorm(1.644854, lower.tail=FALSE),
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, "two.sided"), 
    pnorm(1.644854, lower.tail=FALSE)*2,
    tolerance = 0.05
  )
  expect_equal(
    perm.test(normData, 1.644854, 1, "two.sided"), 
    pnorm(-1.644854, lower.tail=TRUE)*2,
    tolerance = 0.05
  )
})

