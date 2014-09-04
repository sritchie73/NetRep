context("Testing C++ S4 Method Extensions for big.matrix")

# Complete Adjacency
adj <- matrix(rnorm(100), 10, 10)
adjPtr <- as.big.matrix(adj)

# Non square matrix
ns <- matrix(rnorm(30), 3)
nsPtr <- as.big.matrix(ns)
ns2 <- matrix(rnorm(30), 10)
ns2Ptr <- as.big.matrix(ns2)

test_that("scaleBigMatrix is correct", {
  m <- matrix(1:9, 3)
  bigm <- as.big.matrix(m, type="double")
  expect_equivalent(scaleBigMatrix(bigm)[,], scale(m))
})

test_that("checkFinite is correct", {
  set.seed(1)
  options(bigmemory.typecast.warning=FALSE)
  m1 <- matrix(sample(c(1:4, NA), 9, TRUE), 3)
  m2 <- matrix(sample(c(1:4, Inf), 9, TRUE), 3)
  m3 <- matrix(sample(c(1:4, NaN), 9, TRUE), 3)
  m4 <- matrix(sample(c(1:4), 9, TRUE), 3)
  types <- c("char", "short", "integer", "double")
  for (i in 1:4) {
    bigm1 <- as.big.matrix(m1, type=types[i])
    bigm4 <- as.big.matrix(m4, type=types[i])
    expect_error(checkFinite(bigm1))
    expect_null(checkFinite(bigm4))
  }
  bigm2 <- as.big.matrix(m2, type="double")
  bigm3 <- as.big.matrix(m3, type="double")
  expect_error(checkFinite(bigm2))
  expect_error(checkFinite(bigm3))
})