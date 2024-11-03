\name{Rtnmin}
\alias{tn}
\alias{tnbc}
\alias{Rtnmin}
\encoding{UTF-8}
\title{Truncated Newton function minimization}
\description{An R implementation of the truncated Newton method
for minimization of nonlinear functions optionally subject 
to bounds (box) constraints.  The algorithm is based on Nash (1979) 
}
\usage{
   tn(x, fgfun, trace=0, control=list(), \dots)

   tnbc(x, fgfun, lower, upper, trace=0, control=list(), \dots)
}
\arguments{
 \item{x}{A numeric vector of starting estimates.}
 \item{fgfun}{A function that returns the value of the objective at
  the supplied set of parameters \code{par} using auxiliary data in
  \dots. The gradient is returned as attribute "gradient". 
   The first argument of \code{fgfun} must be \code{par}. }
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{trace}{Set > 0 to cause intermediate output to allow progress
       to be followed.}
 \item{control}{A list of control parameters. See the R code for \code{tn()}
     and \code{tnbc()} for details. Most users can omit this.}       
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  Function \code{fgfun} must return a numeric value in list item \code{f}
  and a numeric vector in list item \code{g}.

  These routines can be most easily called using the \code{optimr} 
  function with method "Rtnmin".

  Note that as at 2024-10-23, the evaluation limits maxit and maxfun
  which appear in the code do NOT appear to be active.
}

\value{
  A list with components:
  \item{xstar}{The best set of parameters found.}
  \item{f}{The value of the objective at the best set of parameters found.}
  \item{g}{The gradient of the objective at the best set of parameters found.}
  \item{ierror}{An integer indicating the situation on termination. \code{0}
   indicates that the method believes it has succeeded; \code{2} that
   more than \code{maxfun} (default 150*n, where there are n parameters);
   \code{3} if the line search appears to have failed (which may not be serious);
   and \code{-1} if there appears to be an error in the input parameters.}
  \item{nfngr}{A number giving a measure of how many conjugate gradient solutions
	were used during the minimization process.}
}
\references{

   Stephen G. Nash (1984) "Newton-type minimization via the Lanczos method",
   SIAM J Numerical Analysis, vol. 21, no. 4, pages 770-788.

   For Matlab code, see http://www.netlib.org/opt/tn
}
\seealso{\code{\link{optimr}}}
\examples{
## See tn.Rd

}

\keyword{nonlinear}
\keyword{optimize}

