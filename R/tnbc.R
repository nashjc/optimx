tnbc <- function (x, fgfun, lower, upper, trace=0, ...) {
##---------------------------------------------------------
## this routine solves the optimization problem
##
##   minimize     f(x)
##   subject to   lower <= x <= upper
##
## parameters:
##
## ierror  <-  error code
##             ( 0 => normal return
##             ( 2 => more than maxfun evaluations
##             ( 3 => line search failed (may not be serious)
##             (-1 => error in input parameters
## x        -> initial estimate of the solution; 
## fgfun     -> function routine: in Matlab [f,g] = fgfun(x)
## xstar   <-  the computed solution.
## g       <-  final value of the gradient
## f       <-  final value of the objective function
## lower, upper  -> lower and upper bounds on the variables
##
## this routine sets up all the parameters for lmqnbc:
##
## maxfun - maximum allowable number of function evaluations
## stepmx - maximum allowable step in the linesearch
## accrcy - accuracy of computed function values
## maxit  - maximum number of inner iterations per step
##---------------------------------------------------------
eps    <- .Machine$double.eps
n      <- length(x)
maxit  <- 1 + round((n+1)/2)
maxit  <- min(50, maxit)
maxfun <- 150*n
stepmx <- 10
accrcy <- 100*eps
##---------------------------------------------------------
## return [xstar, f, g, ierror] 
    lmqnbc(x, fgfun, lower, upper, maxit, maxfun, stepmx, accrcy, trace=trace)
}
