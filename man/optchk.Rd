\name{optchk}
\alias{optchk}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  A wrapper function that attempts to check the objective function,
  and optionally the gradient and hessian functions, supplied by the
  user for optimization. It also tries to check the scale of the
  parameters and bounds to see if they are reasonable. 
}
\usage{
optchk(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            control=list(), ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.}
 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  While it can be envisaged that a user would have an analytic hessian but not an analytic
  gradient, we do NOT permit the user to test the hessian in this situation.

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Note that no other attributes of \code{par}
  are copied over. (We have not verified this as at 2009-07-29.)

}
\value{
  A list of the following items:

 \describe{

 \item{grOK}{TRUE if the analytic gradient and a numerical approximation via \code{numDeriv}
   agree within the \code{control$grtesttol} as per the \code{R} code in function 
   \code{grchk}. \code{NULL} if no analytic gradient function is provided. } 
 \item{hessOK}{TRUE if the analytic hessian and a numerical approximation via \code{numDeriv::jacobian}
   agree within the \code{control$hesstesttol} as per the \code{R} code in function 
   \code{hesschk}. NULL if no analytic hessian or no analytic gradient is provided. Note
   that since an analytic gradient must be available for this test, we use the Jacobian of the 
   gradient to compute the Hessian to avoid one level of differencing, though the \code{hesschk}
   function can work without the gradient.}
 \item{scalebad}{TRUE if the larger of the \code{scaleratios} exceeds \code{control$scaletol}}
 \item{scaleratios}{A vector of the parameter and bounds scale ratios. See the function code
   of \code{scalechk} for the computation of these values.}
 }
} 
\references{

 See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

 Nash JC (2014). On Best Practice Optimization Methods in R., 
        \emph{Journal of Statistical Software}, 60(2), 1-14.,
        URL http://www.jstatsoft.org/v60/i02/.

}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
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

myctrl<- ctrldefault(2)
myctrl$trace <- 3
mychk <- optchk(par=c(-1.2,1), fr, grr, lower=rep(-10,2), upper=rep(10,2), control=myctrl)
cat("result of optchk\n")
print(mychk)

}
\keyword{nonlinear}
\keyword{optimize}
