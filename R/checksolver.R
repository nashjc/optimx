checksolver <- function(method, allmeth, allpkg){
#    basestats <- c("Nelder-Mead","BFGS","L-BFGS-B","CG","SANN", "nlm", "nlminb", "hjn")
#     cat("method = ",method,"\n")
#     cat("allmeth:"); print(allmeth)
#     cat("allpkg:");print(allpkg)
#    if (method %in% basestats) return(method)
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
        warning("Package ",pkg," for method ",method," is not available")
        return("????")
      }
    }     
    NULL # just in case
}

checkallsolvers <- function() {
  badmeth <- c() # initially empty
  cc <- ctrldefault(4) # 4 is arbitrary
  ameth <- cc$allmeth
  apkg  <- cc$allpkg
  for (m1 in ameth) {
    p1 <- apkg[which(ameth == m1)]
    cat("Check if method ",m1," from package ",p1," is available\n")
    csres <- checksolver(m1, ameth, apkg)
    if (csres != m1) {
      badmeth <- c(badmeth, m1)  
    }
  }
  if (length(badmeth) > 0) {
    warning("Some methods unavailable -- see badmeth")
  }
  badmeth
}

