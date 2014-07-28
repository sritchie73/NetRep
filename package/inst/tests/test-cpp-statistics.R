context("Testing C++ Implementations of statistics for validity")

# Complete Adjacency
set.seed(1)
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Sparse Adjacency
spr <- adj
spr[sample(1:100, 25)] <- NA
sprPtr <- as.big.matrix(spr)

# Symmetric Matrix
sym <- adj
diag(sym) <- NA
sym[lower.tri(sym)] <- t(sym)[lower.tri(sym)]
symPtr <- as.big.matrix(sym)

test_that("MeanAdj implementation is correct", {
  mod <- sample(1:10, 4)
  expect_equal(meanAdjR(adj, 1:10), meanAdj(adjPtr, 1:10, FALSE))
  expect_equal(meanAdjR(spr, 1:10), meanAdj(sprPtr, 1:10, FALSE))
  expect_equal(meanAdjR(adj, mod), meanAdj(adjPtr, mod, FALSE))
  expect_equal(meanAdjR(sym, mod), meanAdj(symPtr, mod, TRUE))
  expect_error(meanAdj(adjPtr, 0:1, FALSE))
  expect_error(meanAdj(adjPtr, 11:12, FALSE))
})

test_that("kIM implementation is correct", {
  mod <- sample(1:10, 4)
  expect_equal(kIMR(adj, 1:10, TRUE), kIM(adjPtr, 1:10, TRUE))
  expect_equal(kIMR(spr, 1:10, TRUE), kIM(sprPtr, 1:10, TRUE))
  expect_equal(kIMR(adj, mod, TRUE), kIM(adjPtr, mod, TRUE))
  expect_equal(kIMR(adj, mod, FALSE), kIM(adjPtr, mod, FALSE))
})
