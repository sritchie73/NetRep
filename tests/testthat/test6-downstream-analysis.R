context("Testing downstream analysis functions")
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

test_that("'networkProperties' function runs without error", {
  expect_is(
    networkProperties(
      exprSets, coexpSets, adjSets, moduleAssignments, modules="1", verbose=FALSE
    ), "list"
  )
  props <- networkProperties(
    exprSets[[1]][,1:10], coexpSets[[1]][1:10, 1:10], adjSets[[1]][1:10, 1:10],
    verbose=FALSE
  )
  expect_is(props, "list")
})

test_that("'nodeOrder' function runs without error", {
  n <- nodeOrder(
    exprSets, coexpSets, adjSets, moduleAssignments, modules="1", verbose=FALSE
  )
  expect_is(n, "character")
  n <- nodeOrder(
    NULL, coexpSets, adjSets, moduleAssignments, modules=c("1", "7"), 
    orderModules=FALSE, verbose=FALSE
  )
  expect_is(n, "character")
})

test_that("'sampleOrder' function runs without error", {
  s <- sampleOrder(
    exprSets, coexpSets, adjSets, moduleAssignments, modules="1", 
    verbose=FALSE
  )
  expect_is(s, "integer")

  expect_error(
    sampleOrder(
      NULL, coexpSets, adjSets, moduleAssignments, modules=c("1", "7"), 
      simplify=FALSE, verbose=FALSE
    )
  )
})

