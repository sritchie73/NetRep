context("Checking package dependencies")
flag <- FALSE

# First check R version
Rversion <- "3.2.2"
thisVersion <- paste(R.Version()[["major"]], R.Version()[["minor"]], sep=".")
if (thisVersion > Rversion) {
  if (!flag) cat("\n")
  cat(
    "R is newer than NetRep ",
    "(version ", thisVersion, " detected, ",
    "tested with version ", Rversion, ").\n", sep=""
  )
  flag <- TRUE
}

# Check the dependencies
expectedVersionsDepends <- c(
  bigmemory="4.5.8",
  foreach="1.4.3",
  iterators="1.0.8",
  itertools="0.1.3",
  Rcpp="0.12.1",
  utils="3.2.2",
  statmod="1.4.22",
  RcppArmadillo="0.6.200.2.0",
  BH="1.58.0.1",
  RhpcBLASctl="0.15.148",
  abind="1.4.3",
  RColorBrewer="1.1.2"
)

for (pkg in names(expectedVersionsDepends)) {
  pkgVersion <- packageVersion(pkg)
  if (pkgVersion > expectedVersionsDepends[pkg]) {
    if (!flag) cat("\n")
    cat(
      "Package ", pkg, " is newer than NetRep ",
      "(version ", as.character(pkgVersion), " detected, ",
      "tested with version ", expectedVersionsDepends[pkg], ").\n", sep=""
    )
    flag <- TRUE
  }
}

# Check suggested packages
expectedVersionsSuggests <- c(
  testthat="0.11.0",
  doMC="1.3.4",
  doParallel="1.0.10",
  WGCNA="1.48"
)
for (pkg in names(expectedVersionsSuggests)) {
  tryCatch({
    pkgVersion <- packageVersion(pkg)
    if (pkgVersion > expectedVersionsSuggests[pkg]) {
      if (!flag) cat("\n")
      cat(
        "Package ", pkg, " is newer than NetRep ",
        "(version ", as.character(pkgVersion), " detected, ",
        "tested with version ", expectedVersionsSuggests[pkg], ").\n", sep=""
      )
      flag <- TRUE
    }
  }, error = function(e) {
    # ignore
  })
}


# Output message based on test results
if (flag) {
  cat(
    "\nYou can safely ignore these warnings if the following tests pass", 
    "without warnings or errors.\n"
  )
} else {
  cat("All ok!")
}


