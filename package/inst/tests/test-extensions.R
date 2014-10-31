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
  on.exit(unlink(c("*.bin", "*.desc")))
  m <- matrix(1:9, 3)
  as.big.matrix(
    m, type="double", backingfile="tmp.bin", descriptorfile="tmp.desc"
  )
  s <- scaleBigMatrix("tmp.desc", ".")
  scaled <- attach.big.matrix(s)
  expect_equivalent(scaled[,], scale(m))
})

test_that("checkFinite is correct", {
  on.exit(unlink(c("*.bin", "*.desc")))
  set.seed(1)
  options(bigmemory.typecast.warning=FALSE)
  m1 <- matrix(sample(c(1:4, NA), 9, TRUE), 3)
  m2 <- matrix(sample(c(1:4, Inf), 9, TRUE), 3)
  m3 <- matrix(sample(c(1:4, NaN), 9, TRUE), 3)
  m4 <- matrix(sample(c(1:4), 9, TRUE), 3)
  types <- c("char", "short", "integer", "double")
  for (i in 1:4) {
    as.big.matrix(
      m1, type=types[i], backingfile="tmp1.bin", descriptorfile="tmp1.desc"
    )
    as.big.matrix(
      m4, type=types[i], backingfile="tmp4.bin", descriptorfile="tmp4.desc"
    )
    expect_error(checkFinite("tmp1.desc"))
    expect_null(checkFinite("tmp4.desc"))
  }
  bigm2 <- as.big.matrix(
    m2, type="double", backingfile="tmp2.bin", descriptorfile="tmp2.desc"
  )
  bigm3 <- as.big.matrix(
    m3, type="double", backingfile="tmp3.bin", descriptorfile="tmp3.desc"
  )
  expect_error(checkFinite("tmp2.desc"))
  expect_error(checkFinite("tmp3.desc"))
})