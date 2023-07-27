\name{bmstep}
\alias{bmstep}
\encoding{UTF-8}
\title{Compute the maximum step along a search direction.}
\description{Nonlinear optimization problems often have explicit or implicit upper and
lower bounds on the parameters of the function to be miminized or maximized. These are 
called bounds or box constraints. Some of the parameters may be fixed for a given problem
or for a temporary trial. These fixed, or masked, paramters are held at one value during 
a specific 'run' of the optimization.

The bmstep() function computes the maximum step possible (which could be infinite) 
along a particular search direction from current parameters to bounds.

}
\usage{
   bmstep(par, srchdirn, lower=NULL, upper=NULL, bdmsk=NULL, trace=0)
}
\arguments{
 \item{par}{A numeric vector of starting values of the optimization function parameters.}
 \item{srchdirn}{A numeric vector giving the search direction.}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization. Partly for historical reasons, we use the 
     same array during the progress of optimization as an indicator that a 
     parameter is at a lower bound (bdmsk element set to -3) or upper bound (-1).}
 \item{trace}{An integer that controls whether diagnostic information is displayed.
              A positive value displays information, 0 (default) does not.}
}
\details{
     The bmstep function will compute and return (as a double or Inf) the 
     maximum step to the bounds.
}

\value{
  A double precision value or Inf giving the 
     maximum step to the bounds.
}
\examples{
#####################
xx <- c(1, 1)
lo <- c(0, 0)
up <- c(100, 40)
sdir <- c(4,1)
bm <- c(1,1) # both free
ans <- bmstep(xx, sdir, lo, up, bm, trace=1)
# stepsize
print(ans)
# distance
print(ans*sdir)
# New parameters
print(xx+ans*sdir)

}

\keyword{nonlinear}
\keyword{optimize}
\keyword{upper}
\keyword{lower}
\keyword{bound}
\keyword{mask}
