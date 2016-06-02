context("Testing C++ S4 Method Extensions for big.matrix")

test_that("checkFinite is correct", {
  options(bigmemory.typecast.warning=FALSE)
  set.seed(1)
  m1 <- matrix(sample(c(1:2, NA), 9, TRUE), 3)
  m2 <- matrix(sample(c(1:2, Inf), 9, TRUE), 3)
  m3 <- matrix(sample(c(1:2, NaN), 9, TRUE), 3)
  m4 <- matrix(sample(c(1:2), 9, TRUE), 3)
  types <- c("char", "short", "integer", "double")
  for (i in 1:4) {
    bm1 <- as.bigMatrix(m1, file.path(tempdir(), paste0("tmp1-", i)), type=types[i])
    bm4 <- as.bigMatrix(m4, file.path(tempdir(), paste0("tmp4-", i)), type=types[i])
    expect_error(checkFinite(bm1))
    expect_null(checkFinite(bm4))
  }
  bm2 <- as.bigMatrix(m2, file.path(tempdir(), "tmp2"))
  bm3 <- as.bigMatrix(m3, file.path(tempdir(), "tmp3"))
  expect_error(checkFinite(bm2))
  expect_error(checkFinite(bm3))
})

unlink(file.path(tempdir(), c("tmp*", "scaled*")))