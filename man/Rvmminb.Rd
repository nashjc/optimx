\name{Rvmminb}
\alias{Rvmminb}
\encoding{UTF-8}
\title{Variable metric nonlinear function minimization with bounds constraints}
\description{A bounds-constarined R implementation of a variable metric method for minimization
of nonlinear functions subject to bounds (box) constraints and masks 
(fixed parameters). 

   See manual Rvmmin.Rd for more details and examples.
}
\usage{
   Rvmminb(par, fn, gr, lower, upper, bdmsk, control = list(), \dots)
}
\arguments{
 \item{par}{A numeric vector of starting estimates.}
 \item{fn}{A function that returns the value of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. }
 \item{gr}{A function that returns the gradient of the objective at the
   supplied set of parameters \code{par} using auxiliary data in \dots.
   The first argument of \code{fn} must be \code{par}. This function 
   returns the gradient as a numeric vector.

   Note that a gradient function MUST be provided. See the manual for
   \code{Rvmmin}, which is the usual way \code{Rvmminb} is called. The
   user must take responsibility for errors if \code{Rvmminb} is called
   directly. 
 }
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization.}
 \item{control}{
    An optional list of control settings. See the manual Rvmmin.Rd for 
    details. 

 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{
  This routine is intended to be called from \code{Rvmmin}, which will, if
  necessary, supply a gradient approximation. However, some users will want
  to avoid the extra overhead, in which case it is important to provide an
  appropriate and high-accuracy gradient routine.

  Note that bounds checking, if it is carried out, is done by \code{Rvmmin}.

  Functions \code{fn} must return a numeric value.
}

\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of the objective at the best set of parameters found.}
  \item{counts}{A vector of two integers giving the number of function and gradient evaluations.}
  \item{convergence}{An integer indicating the situation on termination of the function. \code{0}
   indicates that the method believes it has succeeded. Other values:
   \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
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
\seealso{\code{\link{optim}}}
\examples{
## See Rvmmin.Rd

}

\keyword{nonlinear}
\keyword{optimize}

