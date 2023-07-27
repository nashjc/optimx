proptimr <- function(opres){
  # to provide compact output of an optimr result
  obname<-deparse(substitute(opres))
  fname<-attr(opres$value, "fname")
  mname<-attr(opres$value, "method")
  cat("Result ",obname,"(",mname," -> ",fname,") calc. min. =", opres$value,
      " at \n")
  n <- length(opres$par)
  if (length(attr(opres$par, "status"))> 0) schar<- attr(opres$par, "status")
  else schar <- rep(" ",n)

  outpar<-""
  for (i in 1:n) cat(opres$par[i],schar[i],"  ")
  cat("\n")
  if (is.null(opres$scounts)){ nfn <- opres$counts[[1]];   ngr <- opres$counts[[2]]; nhe <- NA}
  else {nfn <- opres$scounts[1]; ngr <- opres$scounts[2]; nhe <- opres$scounts[3] }
  cat("After ",nfn," fn evals, and ",ngr," gr evals and ",nhe," hessian evals\n")
  cat("Termination code is ", opres$convergence,":", opres$message,"\n")
  if (! is.null(opres$grad)) {
    cat("Gradient:")
    print(opres$grad)
  }
  if (! is.null(opres$Hess)) {
    cat("Hessian:\n")
    print(opres$Hess)
  }
  cat("\n-------------------------------------------------\n")
}
# -------------- end proptimr ----------------- #
