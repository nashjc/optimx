optchk <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,  
                             control=list(), ...) {

## Should be run whenever we are not sure parameters and function are
## admissible. 
##
## Inputs:  
# par - a vector of initial values for the parameters 
# fn  - A function to be minimized (or maximized)
# gr  - A function to return (as a vector) the gradient 
# hess- A function to return (as a symmetric matrix) the Hessian of the objective 
# lower, upper - Bounds on the variables
# control - A list of control parameters. 
# ... - further arguments to be passed to fn and gr

## Outputs: ?? failed-checks info.

###############################################################################


# Check parameters are in right form
  if (!is.null(dim(par))) stop("Parameters should be a vector, not a matrix!")
  if (! is.vector(par) ) {
	stop("The parameters are NOT in a vector")
  }
  npar<-length(par)
  if (is.null(control)) control <- ctrldefault(npar)
  optchk<-list() # set up output kust of the checks
  # Check parameters in bounds (090601: As yet not dealing with masks ??)
  infeasible<-FALSE
  if (control$trace > 1) cat("Function has ",npar," arguments\n")
  if (! control$have.bounds) { # Don't do the check if we already know there are no bounds
     if (is.null(control$keepinputpar)) {shift2bound <- TRUE }
     else {shift2bound <- ! control$keepinputpar}
     bc <- bmchk(par, lower=lower, upper=upper, bdmsk=rep(1,npar), trace=control$trace, shift2bound)
     if (! bc$admissible) stop("At least one lower bound is > corresponding upper bound")
     if (infeasible && control$dowarn) warning("Parameters may be out of bounds")
     if (control$trace > 0) {
     cat("Parameter relation to bounds\n")
        print(bc$bchar)
     }
     if (bc$parchanged) {
        if (control$trace > 0) cat("parameters have been moved to nearest bounds\n")
        par <- bc$bvec
     }
  }

  # Check if function can be computed
  checkfn <- fnchk(par, fn, trace=control$trace, ...)
  if (checkfn$infeasible) {
     cat("fnchk exit code and msg:",checkfn$excode," ",checkfn$msg,"\n")
     stop("Cannot evaluate function at initial parameters")
  }
  grOK <- NULL
  hessOK <- NULL
  if (! is.null(gr) && ! is.character(gr)){ # check gradient
     if (is.null(control$grtesttol)) stop("optchk: A control$grtesttol is required")
     grOK <- grchk(par, fn, gr, trace=control$trace, testtol=(.Machine$double.eps)^(1/3), ...) 
     if (control$trace > 0) cat("gradient check OK =",grOK,"\n")
     if (! is.null(hess) ) { # if hessian analytic function provided, then check it
       # Note: we only do this if analytic gradient is provided
       if (is.null(control$hesstesttol)) stop("optchk: A control$hesstesttol is required")
       hessOK <- hesschk(par, fn, gr, hess, trace=control$trace, testtol=control$hesstesttol, ...)
       if (control$trace > 0) cat("hessian check OK =",hessOK,"\n")
     }
  } else if (control$trace > 0) cat("Analytic gradient not made available.\n")

# Scaling check  091219
  scalebad <- FALSE
  srat<-scalechk(par, lower, upper,control$dowarn)
  sratv<-c(srat$lpratio, srat$lbratio)
  if (max(sratv,na.rm=TRUE) > control$scaletol) { 
     warnstr<-"Parameters or bounds appear to have different scalings.\n  This can cause poor performance in optimization. \n  It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
     if (control$dowarn) warning(warnstr)
     scalebad <- TRUE
     if (control$trace > 0) {
        cat("Scale check -- log parameter ratio=",srat$lpratio,"  log bounds ratio=",srat$lbratio,"\n")
     }
  }
  optcheck <- list(grOK = grOK, hessOK = hessOK, scalebad = scalebad, scaleratios = sratv)
} ## end of optchk.R

