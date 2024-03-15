\name{snewton}
\alias{snewton}
\alias{snewtm}
\encoding{UTF-8}
\title{Safeguarded Newton methods for function minimization using R functions.}
\description{These versions of the safeguarded Newton solve the Newton equations with
the R function solve(). In \code{snewton} a backtracking line search is used,
while in \code{snewtonm} we rely on a Marquardt stabilization.
}
\usage{
   snewton(par, fn, gr, hess, control = list(trace=0, maxit=500), \dots)

   snewtm(par, fn, gr, hess, bds, control = list(trace=0, maxit=500))
}
\arguments{
 \item{par}{A numeric vector of starting estimates.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{gr}{A function that returns the gradient of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. This function 
   returns the gradient as a numeric vector.}
 \item{hess}{A function to compute the Hessian matrix. This should be provided as a square,
      symmetric matrix.}
 \item{bds}{Result of \code{bmchk()} for the current problem. Contains lower and upper etc.}
 \item{control}{
    An optional list of control settings.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}

\details{
  \code{snewtm} is intended to be called from \code{optimr()}.

  Functions \code{fn} must return a numeric value. \code{gr} must return a vector.
  \code{hess} must return a matrix. 
   The \code{control} argument is a list. See the code for \code{snewton.R} for completeness.
   Some of the values that may be important for users are:
   \describe{
   \item{trace}{Set 0 (default) for no output, > 0 for diagnostic output
      (larger values imply more output).}
   \item{watch}{Set TRUE if the routine is to stop for user input (e.g., Enter)
       after each iteration. Default is FALSE.}
   \item{maxit}{A limit on the number of iterations (default 500 + 2*n where n is
      the number of parameters). This is the maximum number of gradient evaluations 
      allowed.}
   \item{maxfeval}{A limit on the number of function evaluations allowed 
     (default 3000 + 10*n).}
   \item{eps}{ a tolerance used for judging small gradient norm (default = 1e-07).
   	a gradient norm smaller than (1 + abs(fmin))*eps*eps is considered small 
   	enough that a local optimum has been found, where fmin is the current 
   	estimate of the minimal function value. }
   \item{acctol}{To adjust the acceptable point tolerance (default 0.0001) in the test
       ( f <= fmin + gradproj * steplength * acctol ). This test is used to ensure progress
       is made at each iteration. }
   \item{stepdec}{Step reduction factor for backtrack line search (default 0.2)}
   \item{defstep}{Initial stepsize default (default 1)}
   \item{reltest}{Additive shift for equality test (default 100.0)}
  }
  The (unconstrained) solver \code{snewtonmu} proved to be slower than the bounded solver
  called without bounds, so has been withdrawn.

  The \code{snewton} safeguarded Newton uses a simple line search but no linear solution
  stabilization and has demonstrated POOR performance and reliability. NOT recommended.
}

\value{
  A list with components:
  
  \describe{
  \item{par}{The best set of parameters found.}
  \item{value}{The value of the objective at the best set of parameters found.}
  \item{grad}{The value of the gradient at the best set of parameters found. A vector.}
  \item{hessian}{The value of the Hessian at the best set of parameters found. A matrix.}
  \item{counts}{A vector of 4 integers giving number of Newton equation solutions, the number of function
      evaluations, the number of gradient evaluations and the number of hessian evaluations.}
  \item{message}{A message giving some information on the status of the solution.}
  }

}
\references{ 

  Nash, J C (1979, 1990) Compact Numerical Methods for Computers: Linear
     Algebra and Function Minimisation, Bristol: Adam Hilger. Second
     Edition, Bristol: Institute of Physics Publications.

}
\seealso{\code{\link{optim}}}
\examples{
#Rosenbrock banana valley function
f <- function(x){
return(100*(x[2] - x[1]*x[1])^2 + (1-x[1])^2)
}
#gradient
gr <- function(x){
return(c(-400*x[1]*(x[2] - x[1]*x[1]) - 2*(1-x[1]), 200*(x[2] - x[1]*x[1])))
}
#Hessian
h <- function(x) {
a11 <- 2 - 400*x[2] + 1200*x[1]*x[1]; a21 <- -400*x[1]
return(matrix(c(a11, a21, a21, 200), 2, 2))
}

fg <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  val
}
fgh <- function(x){ #function and gradient
  val <- f(x)
  attr(val,"gradient") <- gr(x)
  attr(val,"hessian") <- h(x)
  val
}

x0 <- c(-1.2, 1)

sr <- snewton(x0, fn=f, gr=gr, hess=h, control=list(trace=1))
print(sr)
# Call through optimr to get correct calling sequence, esp. with bounds
srm <- optimr(x0, fn=f, gr=gr, hess=h, control=list(trace=1))
print(srm)

# bounds constrained example

lo <- rep((min(x0)-0.1), 2)
up <- rep((max(x0)+0.1), 2)
# Call through optimr to get correct calling sequence, esp. with bounds
srmb <- optimr(x0, fn=f, gr=gr, hess=h, lower=lo, upper=up, control=list(trace=1))
proptimr(srmb)


#Example 2: Wood function
#
wood.f <- function(x){
  res <- 100*(x[1]^2-x[2])^2+(1-x[1])^2+90*(x[3]^2-x[4])^2+(1-x[3])^2+
    10.1*((1-x[2])^2+(1-x[4])^2)+19.8*(1-x[2])*(1-x[4])
  return(res)
}
#gradient:
wood.g <- function(x){
  g1 <- 400*x[1]^3-400*x[1]*x[2]+2*x[1]-2
  g2 <- -200*x[1]^2+220.2*x[2]+19.8*x[4]-40
  g3 <- 360*x[3]^3-360*x[3]*x[4]+2*x[3]-2
  g4 <- -180*x[3]^2+200.2*x[4]+19.8*x[2]-40
  return(c(g1,g2,g3,g4))
}
#hessian:
wood.h <- function(x){
  h11 <- 1200*x[1]^2-400*x[2]+2;    h12 <- -400*x[1]; h13 <- h14 <- 0
  h22 <- 220.2; h23 <- 0;    h24 <- 19.8
  h33 <- 1080*x[3]^2-360*x[4]+2;    h34 <- -360*x[3]
  h44 <- 200.2
  H <- matrix(c(h11,h12,h13,h14,h12,h22,h23,h24,
                h13,h23,h33,h34,h14,h24,h34,h44),ncol=4)
  return(H)
}
#################################################
w0 <- c(-3, -1, -3, -1)

wd <- snewton(w0, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=1))
print(wd)

# Call through optimr to get correct calling sequence, esp. with bounds
wdm <- optimr(w0, fn=wood.f, gr=wood.g, hess=wood.h, control=list(trace=1))
print(wdm)

}

\keyword{nonlinear}
\keyword{optimize}

