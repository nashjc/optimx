snewton<-function(par, fn, gr, hess, control=list(trace=0, maxit=500),...) {
## Safeguarded Newton minimizer with backtrack line search
##
##Input
##       - par is the initial set of parameter values in vector
##       - fn is the function we wish to minimize
##       - gr is gradient function (required)
##       - hess is hessian function (required)
##       - control is control list (see ctrldefault.R)
##       - ... (dotargs) is exogenous data used in the function fn
##Output (list) -- Note that this does not perfectly match optim() output!!
##       - xs is the parameter vector at the minimum
##       - fv is the fn evaluated at xs
##       - grd is the gradient value (vector)
##       - Hess is the Hessian matrix (note capitalization of this list element name)
##       - niter is the number of interations needed (gradient and Hessian evals).
##       - counts a list of work measures 
##             niter = number of "iterations" i.e., Newton steps
##             nfn   = number of function evaluations
##             ngr   = number of gradient evaluations
##             nhess = number of hessian evaluations
##       - convcode is a number indicating status at termination (0 means OK)
##       - message is a text string returned to indicate status on termination 

npar <- length(par)
nf <- ng <- nh <- niter <- 0 # counters

ctrldefault <- list( ## ?? this should be replaced with call to ctrldefault()
  trace = 0,
  maxit = 500,
  maxfeval = npar*500,
  acctol = 0.0001,
  epstol = .Machine$double.eps,
  stepredn = 0.2, 
  stepmax = 5,
  stepmin = 0,
  offset = 100.0,
  defstep=1,
  bigval = .Machine$double.xmax*0.01,
  watch = FALSE
)  

ncontrol <- names(control)
nctrld <- names(ctrldefault)
for (onename in nctrld) {
  if (! (onename %in% ncontrol)) {
    control[onename]<-ctrldefault[onename]
  }
}
trace <- control$trace # convenience
  out<-NULL # Need this here in case of error return
  if (trace > 2) cat("snewton for problem in ",npar," parameters\n")
  msg <- "Snewton - no msg"
  xb <- par # best so far
  fbest <- fn(xb, ...)
  nf <- nf + 1 
  repeat { # MAIN LOOP
    if (trace > 0) cat(niter," ",nf," ",ng,"  fbest=",fbest,"\n")
    if (trace > 1) print(xb)
    niter <- niter + 1
    grd<-gr(xb,...) # compute gradient
    ng <- ng + 1
    halt <- FALSE # default is keep going
    # tests on too many counts
    if (niter > control$maxit) {
      if (trace > 0) cat("Too many (",niter,") iterations\n")
      halt <- TRUE
      convcode <- 1
      break
    }
    if (nf > control$maxfeval){
      msg <- paste("Too many (",nf," function evaluations")
      if (trace > 0) cat(msg,"\n")
      halt <- TRUE
      convcode <- 91 # value choice?
      break
    }
    gmax <- max(abs(grd))
    if (trace > 1) cat("current gradient norm =",gmax,"\n")
    if (gmax <= control$epstol) {
      msg <- paste("Small gradient norm")
      if (trace > 0) cat(msg,"\n")
      halt <- TRUE
      convcode <- 0 # OK
      break
    }
    # Note if we get here, 
#    if (trace > 0) {cat("Iteration ",niter,":")}
    H<-hess(xb,...)
    nh <- nh + 1
#    eH <- eigen(H, symmetric=TRUE)
#    if (min(eH$values) < 0.1*max(abs(H))*.Machine$double.eps) {       
      d<-try(solve(H, -grd)) # tried change to solve.qr as solve.default seems
          # to give trouble (i.e., error that is not trappable
      if (trace > 3) print(d)
 #   } 
#    else { # failure
#    cat("test inherits try error\n")
    if (inherits(d, "try-error")) {
          msg <- "Failure of default solve of Newton equations"
          if (trace > 0) cat("msg","\n")
          out$par<-xb
          out$value<-fbest
          out$grad<-grd
          out$hessian<-H
          out$counts <- list(niter=niter,  nfn=nf, ngr=ng, nhess=nh)
          out$convcode <- 9030
          out$message <- msg 
#          cat("out"); print(out)         
          return(out)
    }
    if (trace > 2) {
         cat("Search vector:")
         print(d)
    }
    gprj <- as.numeric(crossprod(d, grd))
    if (! is.finite(gprj)) stop("Cannot compute gradient projection")
    if (trace > 0) cat("Gradient projection = ",gprj,"\n")
# Do we want a test on size of gprj?
    st <- control$defstep
    if (gprj > 0) st <- -st # added 180330 to allow downhill
    xnew <- xb + st*d # new point
    if (all((control$offset+xnew) == (control$offset+xb))) {
        convcode <- 92 # no progress
        msg <- "No progress before linesearch!"
        if (trace > 0) cat(msg,"\n")
        break # finished        
    }
    fval <- try(fn(xnew, ...))
    nf <- nf + 1
    if (inherits(fval, "try-error")) stop("snewton: function evaluation error")
    if (trace > 1) {
       cat("f(xnew)=",fval," at ")
       print(xnew)
    }
    if (fval > control$bigval) {
       msg <- "snewton: New function value infinite"
       if (trace > 1) cat(msg,"\n")
#       fval <- control$bigval
       convcode <- 9999
       break
    }
    while ((fval > fbest + control$acctol*abs(st)*gprj) # 180330 Note abs(st)
           && (all((control$offset+xnew) != (control$offset+xb)))) { 
        # continue until satisfied
        st <- st * control$stepredn
        if (trace > 1) cat("Stepsize now =",st,"\n")
        xnew <- xb + st*d # new point
        fval <- fn(xnew, ...)    
        nf <- nf + 1
        if (trace > 1) cat("* f(xnew)=",fval,"\n")
        if (fval > control$bigval) {
           stop("snewton: Function value infinite in backtrack")
        }
    } # end while
    if (all((control$offset+xnew) == (control$offset+xb))) {
        convcode <- 93 # no progress in linesearch
        msg <- "No progress in linesearch!"
        if (trace > 0) cat(msg,"\n")
        break
    }
    if (trace > 1) cat("end major loop\n")  
    xb <- xnew
    fbest <- fval
    if (control$watch) { tmp <- readline("end iteration") }
  } # end repeat
  out$par<-xb
  out$value<-fbest
  out$grad<-grd
  out$hessian<-H
  out$counts <- list(niter=niter,  nfn=nf, ngr=ng, nhess=nh)
  out$convcode <- convcode
  out$message <- msg
  out
}
