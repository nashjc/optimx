\name{ncg}
\alias{ncg}
\encoding{UTF-8}
\title{An R implementation of a  Dai / Yuan nonlinear conjugate gradient algorithm.}
\description{
  Attempts to minimize an unconstrained or bounds (box) and mask constrained function 
  of many parameters by a nonlinear conjugate gradients method using the Dai / Yuan
   update and restart. Based on Nash (1979) Algorithm 22 for its main structure,
   which is method "CG" of the \code{optim()} function that has never performed well.
  Bounds (or box) constraints and masks (equality constraints) can be imposed on 
  parameters. 

  This code is entirely in R to allow users to explore and understand 
  the method. However, \code{ncg()} is intended to be called via \code{optimx::optimr()} and NOT
  called directly, as it has limited sanity checks on the problem provided, since
  such checks are in the \code{optimr()} code.

The earlier \code{Rcgmin()} function
does have such checks, and was originally part of a separate package of the same
name. \code{Rcgmin()} can also be called via \code{optimr()}. It may give slightly
different results due to minor internal coding changes, and is kept available for
backward compatibility with other packages.
}
\usage{
   ncg(par, fn, gr, bds, control = list())
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
   
   Note that \code{ncg} is to be called from \code{optimr} and does 
   NOT allow dot arguments. It is intended to use the internal functions 
   \code{efn} and \code{egr} generated inside \code{optimr()} along with 
   bounds information from \code{bmchk()} available there. 

  The \code{control} argument is a list. See the documentation of \code{ctrldefault()}.

  The source codes \code{Rcgmin} and \code{ncg} for R are still a work 
  in progress, so users should watch the console output. The author welcomes feedback.    
}
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of the objective at the best set of parameters found.}
  \item{counts}{A two-element integer vector giving the number of calls to
          'fn' and 'gr' respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to 'fn'
          to compute a finite-difference approximation to the gradient.}
  \item{convergence}{An integer code. 
	 '0' indicates successful convergence.
         '1' indicates that the function evaluation count 'maxfeval' was reached.
         '2' indicates initial point is infeasible.}
  \item{message}{A character string giving any additional information returned
          by the optimizer, or 'NULL'.}
  \item{bdmsk}{Returned index describing the status of bounds and masks at the
        proposed solution. Parameters for which bdmsk are 1 are unconstrained
        or "free", those with bdmsk 0 are masked i.e., fixed. For historical
        reasons, we indicate a parameter is at a lower bound using -3 
        or upper bound using -1.}
}
\references{

       Dai, Y. H. and Y. Yuan (2001). An efficient hybrid conjugate 
       gradient method for unconstrained optimization. Annals of 
       Operations Research 103 (1-4), 33–47.

       Nash JC (1979). Compact Numerical Methods for Computers: Linear 
       Algebra and Function Minimisation. Adam Hilger, Bristol. Second 
       Edition, 1990, Bristol: Institute of Physics Publications.

       Nash, J. C. and M. Walker-Smith (1987). Nonlinear Parameter 
       Estimation: An Integrated System in BASIC. New York: Marcel Dekker. 
       See https://www.nashinfo.com/nlpe.htm for a downloadable version 
       of this plus some extras.

}
\seealso{\code{\link{optim}}}
\examples{
#####################
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
# Call is from optimr()
ansrosenbrock0 <- optimr(fn=fr,gr=grr, par=c(1,2), method="ncg")
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
onepar<-optimr(startx,sqtst, gr="grfwd", method="ncg", control=list(trace=1)) 
print(onepar)
}

\keyword{nonlinear}
\keyword{optimize}

