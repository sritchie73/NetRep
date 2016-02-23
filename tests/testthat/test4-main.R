context("Testing 'modulePreservation' function")
set.seed(37)
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

# Only analyse modules with > 2 genes
modules <- moduleAssignments[[1]][intersect(names(moduleAssignments[[1]]), 
                                            colnames(adjSets[[2]]))]
modules <- table(modules)
modules <- names(modules[modules > 2])

test_that("Main routine runs and produces sane output", {
  res1 <- modulePreservation(
    exprSets, coexpSets, adjSets, moduleAssignments, modules,
    discovery=1, test=2, nPerm=10, verbose=FALSE, nCores=1
  )
  expect_equal(dim(res1$nulls), c(7, 7, 10))
  expect_equal(dim(res1$observed), c(7, 7))
  expect_equal(dim(res1$p.values), c(7, 7))
  expect_equal(length(res1$propVarsPresent), 7)
  expect_equal(length(res1$nVarsPresent), 7)
  res2 <- modulePreservation(
    NULL, coexpSets, adjSets, moduleAssignments,
    modules, discovery=1, test=2, nPerm=10, 
    verbose=FALSE, nCores=1
  )
  expect_equal(dim(res2$nulls), c(7, 4, 10))
  expect_equal(dim(res2$observed), c(7, 4))
  expect_equal(dim(res2$p.values), c(7, 4))
  expect_equal(length(res2$propVarsPresent), 7)
  expect_equal(length(res2$nVarsPresent), 7)
  
  res1 <- modulePreservation(
    exprSets, coexpSets, adjSets, moduleAssignments, 
    modules, discovery="a", test="b", nPerm=10,
    verbose=FALSE, nCores=1
  )
})
rm(exprSets, coexpSets, adjSets)
gc()
unlink(file.path(tempdir(), 'tmp*'))
