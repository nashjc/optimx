checksolver <- function(method, allmeth, allpkg){
#    Checks if method is available in allmeth
    imeth <- which(method == allmeth)
    if (length(imeth) < 1) {
       warning("Package ",method," not found")
       return(NULL)
    } else {
      pkg <- allpkg[imeth][[1]]
      if ( requireNamespace(pkg, quietly = TRUE)) {
        return(method)
      } else { 
        warning("Package ",pkg," for method ",method,"("," is not available")
        return("????")
      }
    }     
    NULL # just in case
}

checkallsolvers <- function() {
  badmeth <- c() # initially empty
  cc <- ctrldefault(4) # 4 is arbitrary
  ameth <- cc$allmeth
  tnam <- cc$truename
  apkg  <- cc$allpkg
  for (m1 in ameth) {
    imeth <- which(ameth == m1)
    p1 <- apkg[imeth]
    csres <- checksolver(m1, ameth, apkg)
    if (csres != m1) {
      badmeth <- c(badmeth, m1)  
      cat("method ",m1," is missing\n")
    }
    else {
      cat("method ",m1,"(",tnam[imeth],") from package ",p1," is available\n")
    }
  }
  if (length(badmeth) > 0) {
    warning("Some methods unavailable -- see badmeth")
  }
  badmeth
}

