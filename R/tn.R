tn <- function(x, fgfun, trace=0, ...) {
  ## ---------------------------------------------------------
  ##  this routine solves:  minimize f(x)
  ## 
  ##  parameters:
  ## 
  ##  ierror <-  error code
  ##             ( 0 <-> normal return)
  ##             ( 2 <-> more than maxfun evaluations)
  ##             ( 3 <-> line search failed (may not be serious)
  ##             (-1 <-> error in input parameters)
  ##  x       -> initial estimate of the solution; 
  ##  fgfun    -> function routine: [f,g] <- fgfun(x)
  ##  xstar  <-  computed solution.
  ##  g      <-  final value of the gradient
  ##  f      <-  final value of the objective function
  ## 
  ##  This function sets up the parameters for lmqn. They are:
  ##  maxfun - maximum allowable number of function evaluations
  ##  stepmx - maximum allowable step in the linesearch
  ##  accrcy - accuracy of computed function values
  ##  maxit  - maximum number of inner iterations per step
  ## ---------------------------------------------------------
  n  <- length(x)
  maxit <- 1 + round((n+1)/2)
  maxit <- min(50, maxit);
  maxfun <- 150*n;
  stepmx <- 10;
  eps<- .Machine$double.eps
  accrcy <- 100*eps;
  ## result is list of [xstar, f, g, ierror]
  result <- lmqn (x, fgfun, maxit, maxfun, stepmx, accrcy, trace, ...)
  rm("envjn", envir=globalenv())
  result
}


