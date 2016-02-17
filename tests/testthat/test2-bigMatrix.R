context("Testing bigMatrix implementation")
suppressMessages(library(bigmemory))

# Test implementation of `[` with every combination or arguments and matrix types
sampleTF <- function(n) {
  c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)[seq_len(n)]
}
stringifyVector <- function(vec) {
  paste0("c(", paste(vec, collapse=", "), ")")
}

for (mm in 1:3) { # Matrix, single column, and single row behave differently
  switch(mm, 
    "1" = { m1 <- matrix(1:9.0, nrow=3) },
    "2" = { m1 <- matrix(1:9.0, nrow=9) },
    "3" = { m1 <- matrix(1:9.0, ncol=9) })
  for (dd in 1:4) { # no dimnames, rownames only, colnames only, both
    switch(dd, 
      "1" = { dimnames(m1) <- NULL },
      "2" = { dimnames(m1) <- list(letters[seq_len(nrow(m1))], NULL) },
      "3" = { dimnames(m1) <- list(NULL, LETTERS[seq_len(ncol(m1))]) },
      "4" = { dimnames(m1) <- list(letters[seq_len(nrow(m1))], LETTERS[seq_len(ncol(m1))]) })
    save.as.bigMatrix(m1, file.path(tempdir(), "tmp1"), type="integer")
    d1 <- load.bigMatrix(file.path(tempdir(), "tmp1"))
    for (ii in 1:4) {
      # We need to handle cases where there is only 1 row
      switch(ii,
        "1"= {ii <- "1"},
        "2"= { if (mm == 3) { next } else { ii <- "1:2" } },
        "3"= { ii <- ""},
        "4"= { ii <- stringifyVector(sampleTF(nrow(m1))) })
      for (jj in 1:4) {
        # We need to handle cases where there is only 1 column
        switch(jj,
               "1"= { jj <- "1"},
               "2"= { if (mm == 2) { next } else { jj <- "1:2" } },
               "3"= { jj <- ""},
               "4"= { jj <- stringifyVector(sampleTF(ncol(m1))) })
        for (kk in c("", ",drop=FALSE", ",drop=NA", ',drop="a"', ",drop=NULL")) {
          expr1 <- paste0("m1[", ii, ",", jj, kk, "]")
          expr2 <- paste0("d1[", ii, ",", jj, kk, "]")
          # print(paste(expr1, 'vs', expr2, "for mm ==", mm, "dd == ", dd))
          expect_identical(eval(parse(text=expr1)), eval(parse(text=expr2)))
        }
      }
    }
    # Test casting to and from big.matrix
    d2 <- as.big.matrix(d1)
    expect_identical(d2[,], d1[,])
    # Test loading from a "big.matrix"
    d3 <- load.bigMatrix(file.path(tempdir(), "tmp1"))
    expect_identical(d3, d1)
    
    # Check reading and writing
    write.bigMatrix(
      d1, file.path(tempdir(), "tmp1.txt"), 
      sep="\t", quote=FALSE, 
      row.names=ifelse(is.null(rownames(d1)), FALSE, TRUE),
      col.names=ifelse(is.null(colnames(d1)), FALSE, TRUE)
    )
    d5 <- read.bigMatrix(
      file.path(tempdir(), "tmp1.txt"), 
      file.path(tempdir(), "tmp1-1"), 
      type="integer", sep="\t",
      row.names=ifelse(is.null(rownames(d1)), FALSE, TRUE),
      header=ifelse(is.null(colnames(d1)), FALSE, TRUE)
    )
    expect_identical(d1[,], d5[,])
  }
}

cleanTempDir()
