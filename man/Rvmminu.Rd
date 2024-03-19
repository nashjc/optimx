\name{Rvmminu}
\alias{Rvmminu}
\encoding{UTF-8}
\title{Variable metric nonlinear function minimization, unconstrained}
\description{An R implementation of a variable metric method for minimization
of unconstrained nonlinear functions.

See the manual Rvmmin.Rd for details.
}
\usage{
   Rvmminu(par, fn, gr, control = list(), \dots)
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
   \code{Rvmmin}, which is the usual way \code{Rvmminu} is called. The
   user must take responsibility for errors if \code{Rvmminu} is called
   directly. 
 }
 \item{control}{
    An optional list of control settings. See the manual Rvmmin.Rd for 
    details. Some control elements apply only when parameters are bounds
    constrained and are not used in this function.
 }
 \item{\dots}{Further arguments to be passed to \code{fn}.}
}
\details{

  This routine is intended to be called from \code{Rvmmin}, which will, if
  necessary, supply a gradient approximation. However, some users will want
  to avoid the extra overhead, in which case it is important to provide an
  appropriate and high-accuracy gradient routine.

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
}
\seealso{\code{\link{optim}}}
\examples{
####in Rvmmin.Rd ####
}

\keyword{nonlinear}
\keyword{optimize}

