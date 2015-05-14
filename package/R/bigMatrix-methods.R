
setMethod("show", signature(object = "bigMatrix"), function(object) {
  cat(
    'An object of class "bigMatrix"\n',
    'Backingfile:', gsub(".desc", ".bin", object@descriptor), "\n",
    'Current state:', ifelse(object@attached, "attached", "detached"), "\n"
  )
  if(object@attached) {
    cat(" Address:\n")
    show(object@matrix@address)
    cat("\n")
  } 
})