\name{nvm}
\alias{nvm}
\encoding{UTF-8}
\title{Variable metric nonlinear function minimization, driver.}
\description{A driver to call an 
R implementation of a variable metric method for minimization
of nonlinear functions, possibly subject to bounds (box) constraints and masks 
(fixed parameters). The algorithm is based on Nash (1979) Algorithm 21 for main structure,
which is itself drawn from Fletcher's (1970) variable metric code. This is also the basis
of \code{optim()} method 'BFGS' which, however, does not deal with bounds or masks, or 
\code{Rvmmin}. In this method, an approximation to the inverse Hessian (B) is used to generate a 
search direction t = - B \%*\% g, a simple backtracking line search is used until an
acceptable point is found, and the matrix B is updated using a BFGS formula. If no
acceptable point can be found, we reset B to the identity i.e., the search direction
becomes the negative gradient. If the search along the negative gradient is unsuccessful,
the method terminates. 

The code is entirely in R to allow users to explore and understand the method.
However, \code{nvm()} is intended to be called via \code{optimx::optimr()} and NOT
  called directly, as it has limited sanity checks on the problem provided, since
  such checks are in the \code{optimr()} code. 

The earlier \code{Rvmmin()} function
does have such checks, and was originally part of a separate package of the same
name. \code{Rvmmin()} can also be called via \code{optimr()}. It may give slightly
different results due to minor internal coding changes, and is kept available for
backward compatibility with other packages.
}
\usage{
   nvm(par, fn, gr, bds, control = list())
}
\arguments{
 \item{par}{A numeric vector of starting estimates.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par}. This function is created within 
   \code{optimr()} and subsumes auxiliary data in \dots supplied to that wrapper.}
 \item{gr}{A function that returns the gradient of the objective at the
   supplied set of parameters \code{par}. Note that this is usually generated within
   the \code{optimr()} wrapper, where a gradient function or a reference to one of 
   the derivative approximation routines must be provided. See the documentation for
   \code{optimr()} for details.}
 \item{bds}{A list of information resulting from function \code{bmchk} giving
     information on the status of the parameters and bounds and masks.}     
 \item{control}{An optional list of control settings.}
}
\details{
   Function \code{fn} must return a numeric value.

   Note that \code{nvm} is to be called from \code{optimr} and does 
   NOT allow dot arguments. It is intended to use the internal functions 
   \code{efn} and \code{egr} generated inside \code{optimr()} along with 
   bounds information from \code{bmchk()} available there. 

  The \code{control} argument is a list. See the documentation of \code{ctrldefault()}.

  The source codes \code{Rvmmin} and \code{nvm} for R are still a work 
  in progress, so users should watch the console output. The routine
  \code{nvm} attempts to use minimal checking and works only with a
  bounds constrained version of the algorithm, which may work as fast
  as a specific routine for unconstrained problems. This is an open 
  question, and the author welcomes feedback.    
}
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of the objective at the best set of parameters found.}
  \item{counts}{A vector of two integers giving the number of function and gradient evaluations.}
  \item{convergence}{An integer indicating the situation on termination of the function. \code{0}
   indicates that the method believes it has succeeded. Other values:
   \describe{
      \item{\code{0}}{indicates successful termination to an acceptable solution}
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{2}}{indicates that a point with a small gradient norm has been
      found, which is likely a solution.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
   }
  }
  \item{message}{A description of the situation on termination of the function.}
  \item{bdmsk}{Returned index describing the status of bounds and masks at the
        proposed solution. Parameters for which bdmsk are 1 are unconstrained
        or "free", those with bdmsk 0 are masked i.e., fixed. For historical
        reasons, we indicate a parameter is at a lower bound using -3 
         or upper bound using -1.}
}
\references{ 

  Fletcher, R (1970) A New Approach to Variable Metric Algorithms,
     Computer Journal, 13(3), pp. 317-322.

  Nash, J C (1979, 1990) Compact Numerical Methods for Computers: Linear
     Algebra and Function Minimisation, Bristol: Adam Hilger. Second
     Edition, Bristol: Institute of Physics Publications.

}
\seealso{\code{\link{optim}}}
\examples{
# library(optimx)
## Rosenbrock Banana function
fr <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 *      (x2 - x1 * x1))
}
# Call is from optimr(). In this case, small final gradient
ansrosenbrock0 <- optimr(fn=fr,gr=grr, par=c(1,2), method="nvm")
print(ansrosenbrock0) 
#
# Test if 1-parameter problem is possible
#
cat("test n=1 problem using simple squares of parameter\n")

sqtst<-function(xx) {
  res<-sum((xx-2)*(xx-2))
}

nn<-1
startx<-rep(0,nn)
onepar<-optimr(startx,sqtst, gr="grfwd", method="nvm", control=list(trace=1)) 
print(onepar)
}

\keyword{nonlinear}
\keyword{optimize}

