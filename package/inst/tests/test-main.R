context("Testing Main Routine...")
gn1 <- paste0("N_", 1:100)
gn2 <- sample(geneNames1, 100)

coexpSets <- list(
  as.bigMatrix(matrix(rnorm(100*100), 100, dimnames=list(gn1, gn1)), "test1"),
  as.bigMatrix(matrix(rnorm(100*100), 100, dimnames=list(gn2, gn2)), "test2")
)
adjSets <- coexpSets
exprSets <- list(
  as.bigMatrix(matrix(rnorm(50*100), 50, dimnames=list(NULL, gn1)), "dat1"),
  as.bigMatrix(matrix(rnorm(75*100), 75, dimnames=list(NULL, gn2)), "dat2")
)
moduleAssignments <- list(sample(1:7, 100, replace=TRUE), NULL)
names(moduleAssignments[[1]]) <- gn1

# --- 
test_that("Main routine runs and produces sane output", {
  res1 <- netRepMain(
    exprSets, coexpSets, adjSets, moduleAssignments, "1, 2, nPerm=10
  )
  expect_equal(dim(res1$null), c(7, 7, 10))
  expect_equal(dim(res1$observed), c(7, 7))
  expect_equal(dim(res1$p.values), c(7, 7))
  expect_equal(length(res1$propGenesPresent), 7)
  expect_equal(length(res1$genesPresent), 7)
  res2 <- netRepMain(
    exprSets, coexpSets, adjSets, moduleAssignments, 1, 2, nPerm=10
  )
  expect_equal(dim(res2$null), c(7, 4, 10))
  expect_equal(dim(res2$observed), c(7, 4))
  expect_equal(dim(res2$p.values), c(7, 4))
  expect_equal(length(res2$propGenesPresent), 7)
  expect_equal(length(res2$genesPresent), 7)
})

unlink(c("*.bin", "*.desc"))
