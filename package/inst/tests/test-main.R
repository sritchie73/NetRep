context("Testing Main Routine...")

as.big.matrix(matrix(rnorm(100*100), 100, 100),
                   backingfile="test1.bin", descriptorfile="test1.desc")
as.big.matrix(matrix(rnorm(100*100), 100, 100),
                    backingfile="test2.bin", descriptorfile="test2.desc")
as.big.matrix(matrix(rnorm(100*50), 100, 50),
              backingfile="dat1.bin", descriptorfile="dat1.desc")
as.big.matrix(matrix(rnorm(100*75), 100, 75),
              backingfile="dat2.bin", descriptorfile="dat2.desc")
varNameSets <- rep(list(1:100), 2)
nodeNameSets <- rep(list(1:100), 2)
adjSets <- list("test1.desc", "test2.desc")
nodeLabelSets <- list(sample(1:7, 100, replace=TRUE), NULL)
datSets = list("dat1.desc", "dat2.desc")
names(nodeLabelSets[[1]]) <- 1:100

test_that("Main routine runs and produces sane output", {
  res1 <- netRepMain(datSets, varNameSets, adjSets, nodeNameSets,
                     nodeLabelSets, 1, 2, nPerm=10, tailApprox=F)
  expect_equal(dim(res1$null), c(7, 9, 10))
  expect_equal(dim(res1$observed), c(7, 9))
  expect_equal(dim(res1$p.values), c(7, 9))
  expect_equal(length(res1$overlapProp), 7)
  expect_equal(length(res1$overlapSize), 7)
  res2 <- netRepMain(NULL, NULL, adjSets, nodeNameSets,
                     nodeLabelSets, 1, 2, nPerm=10, tailApprox=F)
  expect_equal(dim(res2$null), c(7, 6, 10))
  expect_equal(dim(res2$observed), c(7, 6))
  expect_equal(dim(res2$p.values), c(7, 6))
  expect_equal(length(res2$overlapProp), 7)
  expect_equal(length(res2$overlapSize), 7)
})

unlink(c("*.bin", "*.desc"))
