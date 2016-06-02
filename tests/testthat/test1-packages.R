context("Checking package dependencies")
flag <- FALSE

# First check R version
Rversion <- "3.2.4"
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
testedVersionsDepends <- c(
  foreach="1.4.3",
  Rcpp="0.12.4",
  utils="3.2.4",
  statmod="1.4.24",
  RcppArmadillo="0.6.700.6.0",
  BH="1.60.0.1",
  RhpcBLASctl="0.15.148",
  abind="1.4.3",
  RColorBrewer="1.1.2"
)

for (pkg in names(testedVersionsDepends)) {
  pkgVersion <- packageVersion(pkg)
  if (pkgVersion > testedVersionsDepends[pkg]) {
    if (!flag) cat("\n")
    cat(
      "Package ", pkg, " is newer than NetRep ",
      "(version ", as.character(pkgVersion), " detected, ",
      "tested with version ", testedVersionsDepends[pkg], ").\n", sep=""
    )
    flag <- TRUE
  }
}

# Check suggested packages
testedVersionsSuggests <- c(
  testthat="0.11.0"
)
for (pkg in names(testedVersionsSuggests)) {
  tryCatch({
    pkgVersion <- packageVersion(pkg)
    if (pkgVersion > testedVersionsSuggests[pkg]) {
      if (!flag) cat("\n")
      cat(
        "Package ", pkg, " is newer than NetRep ",
        "(version ", as.character(pkgVersion), " detected, ",
        "tested with version ", testedVersionsSuggests[pkg], ").\n", sep=""
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


