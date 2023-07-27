# hjn.R -- R implementation of J Nash BASIC HJG.BAS 20160705
hjn <- function(par, fn, lower=-Inf, upper=Inf, bdmsk=NULL, control=list(trace=0), ...){
  n <- length(par) # number of parameters
  if (! is.null(control$maximize) && control$maximize) 
           stop("Do NOT try to maximize with hjn()")
  if (is.null(control$trace)) control$trace <- 0 # just in case
  if (is.null(control$stepsize)) {
     stepsize <- 1 # initial step size (could put in control())
  } else { stepsize <- control$stepsize }
  # Want stepsize positive or bounds get messed up
  if (is.null(control$stepredn)) {
     stepredn <- .1 # defined step reduction (again in control()?)
  } else { stepredn <- control$stepredn }
  if (is.null(control$maxfeval)) control$maxfeval<-2000*n
  if (is.null(control$eps)) control$epsl<-1e-07
  steptol <- control$eps
  # Hooke and Jeeves with bounds and masks
  if (length(upper) == 1) upper <- rep(upper, n)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (is.null(bdmsk)) { 
      bdmsk <- rep(1,n)
      idx <- 1:n 
  } else { idx <- which(bdmsk != 0) } # define masks
  if (any(lower >= upper)){
      warning("hjn: lower >= upper for some parameters -- set masks")
      bdmsk[which(lower >= upper)] <- 0
      idx <- which(bdmsk != 0)
  }
  if (control$trace > 0) {
    cat("hjn:bdmsk:")
    print(bdmsk)
#  cat("idx:")
#  print(idx)
  }
  nac <- length(idx)
  if (is.null(control$reltest)) { reltest = 100 } else {reltest <- control$reltest}
  if (any(par < lower) || any(par > upper)) stop("hjn: initial parameters out of bounds")
  pbase <- par # base parameter set (fold is function value)
  f <- fn(par, ...) # fn at base point
  fmin <- fold <- f # "best" function so far
  pbest <- par # Not really needed 
  fcount <- 1 # count function evaluations, compare with maxfeval
#    cat(fcount, "  f=",fold," at ")
#    print(par)
#    tmp <- readline("cont.")
  keepgoing <- TRUE
  ccode <- 1 # start assuming won't get to solution before feval limit
  while (keepgoing) {    
    # exploratory search -- fold stays same for whole set of axes
    if (control$trace > 0) cat("Exploratory move - stepsize = ",stepsize,"\n")
    if (control$trace > 1) {
       cat("p-start:")
       print(par)
    }
    for (jj in 1:nac){ # possibly could do this better in R
       # use unmasked parameters
       j <- idx[jj]
       ptmp <- par[j]
       doneg <- TRUE # assume we will do negative step
       if (ptmp + reltest < upper[j] + reltest) { # Not on upper bound so do pos step 
          par[j] <- min(ptmp+stepsize, upper[j])
          if ((par[j] + reltest) != (ptmp + reltest)) {
             fcount <- fcount + 1
             f <- fn(par, ...)
#               cat(fcount, "  f=",f," at ")
#               print(par)
             if (f < fmin) {
                fmin <- f
                pbest <- par
#                  cat("*")
                doneg <- FALSE # only case where we don't do neg
                resetpar <- FALSE
             } 
#             tmp <- readline("cont>")
          } 
       } # end not on upper bound
       if (fcount >= control$maxfeval) break
       if (doneg) {
         resetpar <- TRUE # going to reset the paramter unless we improve
         if ((ptmp + reltest) > (lower[j] + reltest)) { # can do negative step
            par[j] <- max((ptmp - stepsize), lower[j])
            if ((par[j] + reltest) != (ptmp + reltest)) {
               fcount <- fcount + 1
               f <- fn(par, ...)
#                 cat(fcount, "  f=",f," at ")
#                 print(par)
               if (f < fmin) {
                  fmin <- f
                  pbest <- par
#                  cat("*")
                  resetpar <- FALSE # don't reset parameter
               } 
#              tmp <- readline("cont<")
            }
         } #  neg step possible
       } # end doneg
       if (resetpar) { par[j] <- ptmp }
    } # end loop on axes
    if (fcount >= control$maxfeval) {
        ccode <- 1
        if (control$trace > 0) cat("Function count limit exceeded\n")
        ans <- list(par=pbest, value=fmin, counts=c(fcount, NA), convergence=ccode)
        return(ans)
    }
    if (control$trace > 0) { 
       cat("axial search with stepsize =",stepsize,"  fn value = ",fmin,"  after ",fcount,"  maxfeval =", control$maxfeval,"\n")
    }
    if (fmin < fold) { # success -- do pattern move (control$trace > 0) cat("Pattern move \n")
       if (control$trace > 1) {
          cat("PM from:")
          print(par)
          cat("pbest:")
          print(pbest)
       }
       for (jj in 1:nac) { # Note par is best set of parameters
          j <- idx[jj]
          ptmp <- 2.0*par[j] - pbase[j]
          if (ptmp > upper[j]) ptmp <- upper[j]
          if (ptmp < lower[j]) ptmp <- lower[j]
          pbase[j] <- par[j]
          par[j] <- ptmp 
       }
       fold <- fmin
       if (control$trace > 1) {
          cat("PM to:")
          print(par)
       }
    # Addition to HJ -- test new  base
#       fcount <- fcount + 1
#       f <- fn(par, ...)
#         cat(fcount, "  f=",f," at ")
#         print(par)
#         tmp <- readline("PM point")
#       if (f < fmin) {
#         if (control$trace > 0) {cat("Use PM point as new base\n")}
#         pbest <- pbase <- par
#       }
    } else { # no success in Axial Seart, so reduce stepsize
       if (fcount >= control$maxfeval) {
        ccode <- 1
        if (control$trace > 0) cat("Function count limit exceeded\n")
        ans <- list(par=pbest, value=fmin, counts=c(fcount, NA), convergence=ccode)
        return(ans)
       }
       # first check if point changed
       samepoint <- identical((par + reltest),(pbase + reltest))
       if (samepoint) { 
          stepsize <- stepsize*stepredn
          if (control$trace > 1) cat("Reducing step to ",stepsize,"\n")
          if (stepsize <= steptol) keepgoing <- FALSE
          ccode <- 0 # successful convergence
       } else { # return to old base point
          if (control$trace > 1) {
             cat("return to base at:")
             print(pbase)
          }
          par <- pbase
       }
    }
    if (fcount >= control$maxfeval) {
        ccode <- 1
        if (control$trace > 0) cat("Function count limit exceeded\n")
        ans <- list(par=pbest, value=fmin, counts=c(fcount, NA), convergence=ccode)
        return(ans)
    }
  } # end keepgoing loop 
  if ( control$trace > 1 ) {
    if (identical(pbest, pbase)) {cat("pbase = pbest\n") }
    else { cat("BAD!: pbase != pbest\n") } 
  }
   
  ans <- list(par=pbest, value=fmin, counts=c(fcount, NA), convergence=ccode)
}
