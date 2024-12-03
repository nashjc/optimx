proptimr <- function(opres, nlim=5){
  # to provide compact output of an optimr result
  # 20241202: limit number of output values to nlim, default 5
  obname<-deparse(substitute(opres))
  fname<-attr(opres$value, "fname")
  mname<-attr(opres$value, "method")
  cat("Result ",obname,"(",mname," -> ",fname,") calc. min. =", opres$value,
      " at \n")
  n <- length(opres$par)
  if (length(attr(opres$par, "status"))> 0) schar<- attr(opres$par, "status")
  else schar <- rep(" ",n)

  outpar<-""
  for (i in 1:nlim) cat(opres$par[i],schar[i],"  ")
  cat("\n")
  # Note we use "counts" if "scounts" not available
  if (is.null(opres$scounts)){ nfn <- opres$counts[[1]];   ngr <- opres$counts[[2]]; nhe <- NA}
  else {nfn <- opres$scounts[1]; ngr <- opres$scounts[2]; nhe <- opres$scounts[3] }
  cat("After ",nfn," fn evals, and ",ngr," gr evals and ",nhe," hessian evals\n")
  cat("Termination code is ", opres$convergence,":", opres$message,"\n")
  if (! is.null(opres$grad)) {
    cat("Gradient:")
    print(opres$grad[1:nlim])
  }
  if (! is.null(opres$Hess)) {
    cat("Hessian (first ",nlim," by ",nlim,":\n")
    print(opres$Hess[1:nlim, 1:nlim])
  }
  cat("\n-------------------------------------------------\n")
}
# -------------- end proptimr ----------------- #
