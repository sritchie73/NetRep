context("Testing bigMatrix implementation")
library(bigmemory)

test_that("Testing with no column and rownames", {
  m1 <- matrix(1:9.0, nrow=3)
  save.as.bigMatrix(m1, file.path(tempdir(), "tmp1"), type="integer")
  d1 <- load.bigMatrix(file.path(tempdir(), "tmp1"))
  
  # Test various permutations and combinations of arguments for subsetting
  for (ii in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)")) {
    for (jj in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)")) {
      for (kk in c("", ",drop=FALSE", ",drop=NA", ',drop="a"', ",drop=NULL")) {
        expect_identical(
          eval(parse(text=paste0("m1[", ii, ",", jj, kk, "]"))),
          eval(parse(text=paste0("d1[", ii, ",", jj, kk, "]"))),
        )
      }
    }
  }
  
  # Test casting to and from big.matrix
  d2 <- as.big.matrix(d1)
  expect_identical(d2[,], d1[,])
  # Test loading from a "big.matrix"
  d3 <- load.bigMatrix(file.path(tempdir(), "tmp1"))
  expect_identical(d3, d1)
  # Make sure the row and column names are put in the right files
  d4 <- load.bigMatrix(file.path(tempdir(), "tmp1"))
  expect_identical(d1, d4)
  
  # Check reading and writing
  write.bigMatrix(
    d1, file.path(tempdir(), "tmp1.txt"), 
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
  )
  d5 <- read.bigMatrix(
    file.path(tempdir(), "tmp1.txt"), 
    file.path(tempdir(), "tmp1"), 
    type="integer", sep="\t", 
    header=FALSE, row.names=FALSE
  )
  expect_identical(d1, d5)
})

test_that("Testing with rownames only", {
  m1 <- matrix(1:9, nrow=3, dimnames=list(letters[1:3], NULL))
  save.as.bigMatrix(m1, file.path(tempdir(), "tmp2"), type="integer")
  d1 <- load.bigMatrix(file.path(tempdir(), "tmp2"))

  # Test various permutations and combinations of arguments for subsetting
  for (ii in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)", '"a"', 'c("a","b")')) {
    for (jj in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)")) {
      for (kk in c("", ",drop=FALSE", ",drop=NA", ',drop="a"', ",drop=NULL")) {
        expect_identical(
          eval(parse(text=paste0("m1[", ii, ",", jj, kk, "]"))),
          eval(parse(text=paste0("d1[", ii, ",", jj, kk, "]"))),
        )
      }
    }
  }
  
  # Test casting to and from big.matrix
  d2 <- as.big.matrix(d1)
  expect_identical(d2[,], d1[,])
  # Test loading from a "big.matrix"
  expect_warning(d3 <- load.bigMatrix(file.path(tempdir(), "tmp2"))) 
  expect_identical(d3, d1)
  # Make sure the row and column names are put in the right files
  d4 <- load.bigMatrix(file.path(tempdir(), "tmp2"))
  expect_identical(d1, d4)
  # Check reading and writing
  write.bigMatrix(
    d1, file.path(tempdir(), "tmp2.txt"), 
    sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE
  )
  d5 <- read.bigMatrix(
    file.path(tempdir(), "tmp2.txt"), file.path(tempdir(), "tmp2"), 
    type="integer", sep="\t", header=FALSE, row.names=TRUE
  )
  expect_identical(d1, d5)
})

test_that("Testing with colnames only", {
  m1 <- matrix(1:9, nrow=3, dimnames=list(NULL, letters[1:3]))
  save.as.bigMatrix(m1, file.path(tempdir(), "tmp3"), type="integer")
  d1 <- load.bigMatrix(file.path(tempdir(), "tmp3"))

  # Test various permutations and combinations of arguments for subsetting
  for (ii in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)")) {
    for (jj in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)", '"a"', 'c("a","b")')) {
      for (kk in c("", ",drop=FALSE", ",drop=NA", ',drop="a"', ",drop=NULL")) {
        expect_identical(
          eval(parse(text=paste0("m1[", ii, ",", jj, kk, "]"))),
          eval(parse(text=paste0("d1[", ii, ",", jj, kk, "]"))),
        )
      }
    }
  }
  
  # Test casting to and from big.matrix
  d2 <- as.big.matrix(d1)
  expect_identical(d2[,], d1[,])
  # Test loading from a "big.matrix"
  expect_warning(d3 <- load.bigMatrix(file.path(tempdir(), "tmp3"))) 
  expect_identical(d3, d1)
  # Make sure the row and column names are put in the right files
  d4 <- load.bigMatrix(file.path(tempdir(), "tmp3"))
  expect_identical(d1, d4)
  # Check reading and writing
  write.bigMatrix(
    d1, file.path(tempdir(), "tmp3.txt"), 
    sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE
  )
  d5 <- read.bigMatrix(
    file.path(tempdir(), "tmp3.txt"), file.path(tempdir(), "tmp3"), 
    type="integer", sep="\t", header=TRUE, row.names=FALSE
  )
  expect_identical(d1, d5)
})

test_that("Testing with both row and column names", {
  m1 <- matrix(1:9, nrow=3, dimnames=list(LETTERS[1:3], letters[1:3]))
  save.as.bigMatrix(m1, file.path(tempdir(), "tmp4"), type="integer")
  d1 <- load.bigMatrix(file.path(tempdir(), "tmp4"))
  expect_identical(m1, d1[,])
  
  # Test various permutations and combinations of arguments for subsetting
  for (ii in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)", '"A"', 'c("A","B")')) {
    for (jj in c("1", "1:2", "", "c(TRUE, FALSE, TRUE)", '"a"', 'c("a","b")')) {
      for (kk in c("", ",drop=FALSE", ",drop=NA", ',drop="a"', ",drop=NULL")) {
        expect_identical(
          eval(parse(text=paste0("m1[", ii, ",", jj, kk, "]"))),
          eval(parse(text=paste0("d1[", ii, ",", jj, kk, "]"))),
        )
      }
    }
  }
  
  # Test casting to and from big.matrix
  d2 <- as.big.matrix(d1)
  expect_identical(d2[,], d1[,])
  # Test loading from a "big.matrix"
  expect_warning(d3 <- load.bigMatrix(file.path(tempdir(), "tmp4"))) 
  expect_identical(d3, d1)
  # Make sure the row and column names are put in the right files
  d4 <- load.bigMatrix(file.path(tempdir(), "tmp4"))
  expect_identical(d1, d4)
  write.bigMatrix(
    d1, file.path(tempdir(), "tmp4.txt"), 
    sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE
  )
  d5 <- read.bigMatrix(
    file.path(tempdir(), "tmp4.txt"), file.path(tempdir(), "tmp4"), 
    type="integer", sep="\t", header=TRUE, row.names=TRUE
  )
  expect_identical(d1, d5)
})

unlink(file.path(tempdir(), "tmp*"))
