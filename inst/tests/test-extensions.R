context("Testing C++ S4 Method Extensions for big.matrix")

test_that("scaleBigMatrix is correct", {
  m <- matrix(1:9, 3)
  bm <- as.bigMatrix(m, "tmp1")
  s <- scaleBigMatrix(bm, ".")
  expect_equivalent(s[,], scale(m))
})

test_that("checkFinite is correct", {
  options(bigmemory.typecast.warning=FALSE)
  set.seed(1)
  m1 <- matrix(sample(c(1:2, NA), 9, TRUE), 3)
  m2 <- matrix(sample(c(1:2, Inf), 9, TRUE), 3)
  m3 <- matrix(sample(c(1:2, NaN), 9, TRUE), 3)
  m4 <- matrix(sample(c(1:2), 9, TRUE), 3)
  types <- c("char", "short", "integer", "double")
  for (i in 1:4) {
    bm1 <- as.bigMatrix(m1, "tmp1", type=types[i])
    bm4 <- as.bigMatrix(m4, "tmp4", type=types[i])
    expect_error(checkFinite(bm1))
    expect_null(checkFinite(bm4))
  }
  bm2 <- as.bigMatrix(m2, "tmp2")
  bm3 <- as.bigMatrix(m3, "tmp3")
  expect_error(checkFinite(bm2))
  expect_error(checkFinite(bm3))
})

test_that("rangeBigMatrix is correct", {
  m1 <- matrix(rnorm(30), ncol=10)
  bm1 <- as.bigMatrix(m1, "tmp")
  expect_equivalent(rangeBigMatrix(bm1), range(m1))
  expect_equivalent(rangeBigMatrix(bm1, c(2,4,5,1)), range(m1[,c(2,4,5,1)]))
})

unlink(c("tmp*", "scaled*"))