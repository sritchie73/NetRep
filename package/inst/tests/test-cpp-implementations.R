context("Testing C++ implementations for validity")

# Complete Adjacency
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Sparse Adjacency
spr <- adj
spr[sample(1:100, 25)] <- NA
sprPtr <- as.big.matrix(spr)

test_that("diag and diag<- are correct") {
  expect_that(diag(adjPtr) == diag(adj))
  expect_that(diag(sprPtr) == diag(spr))
  diag(adj) <- 2
  diag(adjPtr) <- 2
  expect_that(all(adjPtr[] == adj))
  diag(adj) <- 1:10
  diag(adjPtr) <- 1:10
  expect_that(all(adjPtr[] == adj))
}

test_that("MeanAdj implementation is correct", {
  mod <- sample(1:10, 4)
  expect_equal(meanAdjR(adj, 1:10), meanAdj(adjPtr, 1:10))
  expect_equal(meanAdjR(spr, 1:10), meanAdj(sprPtr, 1:10))
  expect_equal(meanAdjR(adj, mod), meanAdj(adjPtr, mod))
  expect_error(meanAdj(adjPtr, 0:1))
  expect_error(meanAdj(adjPtr, 11:12))
})

test_that("kIM implementation is correct", {
  mod <- sample(1:10, 4)
  expect_equal(kIMR(adj, 1:10, TRUE), kIM(adjPtr, 1:10, TRUE))
  expect_equal(kIMR(spr, 1:10, TRUE), kIM(sprPtr, 1:10, TRUE))
  expect_equal(kIMR(adj, mod, TRUE), kIM(adjPtr, mod, TRUE))
  expect_equal(kIMR(adj, mod, FALSE), kIM(adjPtr, mod, FALSE))
})