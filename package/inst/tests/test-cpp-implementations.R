context("Testing C++ implementations of statistics for validity")

# Complete Adjacency
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Sparse Adjacency
spr <- adj
spr[sample(1:100, 25)] <- NA
sprPtr <- as.big.matrix(spr)

test_that("MeanAdj implementation is correct", {
  mod <- sample(1:10, 4)
  expect_equal(meanAdjR(adj, 1:10, TRUE), meanAdj(adjPtr, 1:10, TRUE))
  expect_equal(meanAdjR(spr, 1:10, TRUE), meanAdj(sprPtr, 1:10, TRUE))
  expect_equal(meanAdjR(adj, 1:10, FALSE), meanAdj(adjPtr, 1:10, FALSE))
  expect_equal(meanAdjR(adj, mod, TRUE), meanAdj(adjPtr, mod, TRUE))
})