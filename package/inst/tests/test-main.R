context("Testing Main Routine...")
library(netrep)
gn1 <- paste0("N_", 1:100)
gn2 <- paste0("N_", seq(2, 200, length=100))

coexpSets <- list(
  as.bigMatrix(matrix(rnorm(100*100), 100, dimnames=list(gn1, gn1)), "tmp1"),
  as.bigMatrix(matrix(rnorm(100*100), 100, dimnames=list(gn2, gn2)), "tmp2")
)
adjSets <- coexpSets
exprSets <- list(
  as.bigMatrix(matrix(rnorm(50*100), 50, dimnames=list(NULL, gn1)), "tmp3"),
  as.bigMatrix(matrix(rnorm(75*100), 75, dimnames=list(NULL, gn2)), "tmp4")
)
moduleAssignments <- list(sample(1:7, 100, replace=TRUE), NULL)
names(moduleAssignments[[1]]) <- gn1
test_that("Network properties function runs without error", {
  expect_is(
    networkProperties(
      exprSets, coexpSets, adjSets, moduleAssignments, modules="1"
    ), "list"
  )
  expect_is(
    networkProperties(
      exprSets[[1]][,1:10], coexpSets[[1]][1:10, 1:10], adjSets[[1]][1:10, 1:10]
    ), "list"
  )
})

test_that("Main routine runs and produces sane output", {
  res1 <- modulePreservation(
    exprSets, coexpSets, adjSets, moduleAssignments, 1, 2, nPerm=10
  )
  expect_equal(dim(res1$nulls), c(7, 7, 10))
  expect_equal(dim(res1$observed), c(7, 7))
  expect_equal(dim(res1$p.values), c(7, 7))
  expect_equal(length(res1$propGenesPresent), 7)
  expect_equal(length(res1$nGenesPresent), 7)
  res2 <- netRepMain(
    NULL, coexpSets, adjSets, moduleAssignments, 1, 2, nPerm=10
  )
  expect_equal(dim(res2$nulls), c(7, 4, 10))
  expect_equal(dim(res2$observed), c(7, 4))
  expect_equal(dim(res2$p.values), c(7, 4))
  expect_equal(length(res2$propGenesPresent), 7)
  expect_equal(length(res2$nGenesPresent), 7)
})
unlink('tmp*')
