context("Testing C++ implementations of statistics for validity")

# Complete Adjacency
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Sparse Adjacency
spr <- adj
spr[sample(1:100, 25)] <- NA
sprPtr <- as.big.matrix(spr)

test_that("MeanAdj implementation is correct:", {
  mean.diag <- function(x, n) {
    (sum(x) - sum(diag(x))) / (n*n - n)
  }
  expect_equal(mean(adj), FastModPres:::meanAdj(adjPtr, 1:10, TRUE))
  expect_equal(mean(spr, na.rm=TRUE), FastModPres:::meanAdj(sprPtr, 1:10, TRUE))
  expect_equal(mean.diag(adj, 10), FastModPres:::meanAdj(adjPtr, 1:10, FALSE))
})