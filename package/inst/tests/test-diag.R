context("Testing C++ S4 Method Extensions for big.matrix")

# Complete Adjacency
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Non square matrix
ns <- matrix(rnorm(30), 3)
nsPtr <- as.big.matrix(ns)
ns2 <- matrix(rnorm(30), 10)
ns2Ptr <- as.big.matrix(ns2)

test_that("diag is correct", {
  expect_equal(diag(adjPtr), diag(adj))
  expect_equal(diag(nsPtr), diag(nsPtr))
  expect_equal(diag(ns2Ptr), diag(ns2Ptr))
})

test_that("diag<- is correct", {  
  diag(adj) <- 2
  diag(adjPtr) <- 2
  diag(ns) <- 2
  diag(nsPtr) <- 2
  diag(ns2) <- 2
  diag(ns2Ptr) <- 2
  expect_equal(adjPtr[], adj)
  expect_equal(nsPtr[], ns)
  expect_equal(ns2Ptr[], ns2)
  diag(adj) <- 1:10
  diag(adjPtr) <- 1:10
  expect_equal(adjPtr[], adj)
  diag(ns) <- 1:3
  diag(nsPtr) <- 1:3
  expect_equal(nsPtr[], ns)
  diag(ns2) <- 1:3
  diag(ns2Ptr) <- 1:3
  expect_equal(ns2Ptr[], ns2)
})