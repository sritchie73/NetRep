context("Testing 'modulePreservation' function")
gn1 <- paste0("N_", 1:100)
gn2 <- paste0("N_", seq(2, 200, length=100))

coexpSets <- list(
  a=as.bigMatrix(
    matrix(rnorm(100*100), 100, dimnames=list(gn1, gn1)), 
    file.path(tempdir(), "tmp1")
  ),
  b=as.bigMatrix(
    matrix(rnorm(100*100), 100, dimnames=list(gn2, gn2)), 
    file.path(tempdir(), "tmp2")
  )
)
adjSets <- coexpSets
exprSets <- list(
  a=as.bigMatrix(
    matrix(rnorm(50*100), 50, dimnames=list(NULL, gn1)),
    file.path(tempdir(), "tmp3")
  ),
  b=as.bigMatrix(
    matrix(rnorm(75*100), 75, dimnames=list(NULL, gn2)),
    file.path(tempdir(), "tmp4")
  )
)
moduleAssignments <- list(a=sample(1:7, 100, replace=TRUE), b=NULL)
names(moduleAssignments[[1]]) <- gn1

test_that("Main routine runs and produces sane output", {
  res1 <- modulePreservation(
    exprSets, coexpSets, adjSets, moduleAssignments, 1, 2, nPerm=10, 
    keepNulls=TRUE, verbose=FALSE
  )
  expect_equal(dim(res1$nulls), c(7, 7, 10))
  expect_equal(dim(res1$observed), c(7, 7))
  expect_equal(dim(res1$p.values), c(7, 7))
  expect_equal(length(res1$propVarsPresent), 7)
  expect_equal(length(res1$nVarsPresent), 7)
  res2 <- modulePreservation(
    NULL, coexpSets, adjSets, moduleAssignments, 1, 2, nPerm=10,
    keepNulls=TRUE, verbose=FALSE
  )
  expect_equal(dim(res2$nulls), c(7, 4, 10))
  expect_equal(dim(res2$observed), c(7, 4))
  expect_equal(dim(res2$p.values), c(7, 4))
  expect_equal(length(res2$propVarsPresent), 7)
  expect_equal(length(res2$nVarsPresent), 7)
  
  moduleAssignments[[2]]<- sample(1:9, 100, replace=TRUE)
  names(moduleAssignments[[2]]) <- gn2
  res1 <- modulePreservation(
    exprSets, coexpSets, rev(adjSets), moduleAssignments, "a", "b", nPerm=10,
    include=c(1,2,3), keepNulls=TRUE, verbose=FALSE
  )
})
unlink(file.path(tempdir(), 'tmp*'))
