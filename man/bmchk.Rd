\name{bmchk}
\alias{bmchk}
\encoding{UTF-8}
\title{Check bounds and masks for parameter constraints used in nonlinear optimization}
\description{Nonlinear optimization problems often have explicit or implicit upper and
lower bounds on the parameters of the function to be miminized or maximized. These are 
called bounds or box constraints. Some of the parameters may be fixed for a given problem
or for a temporary trial. These fixed, or masked, paramters are held at one value during 
a specific 'run' of the optimization.

It is possible that the bounds are inadmissible, that is, that at least one lower bound
exceeds an upper bound. In this case we set the flag \code{admissible} to FALSE.

Parameters that are outside the bounds are moved to the nearest bound and the flag
\code{parchanged} is set TRUE. However, we DO NOT change masked parameters, and they
may be outside the bounds. This is an implementation choice, since it may be useful
to test objective functions at point outside the bounds.

The package bmchk is essentially a test of the R function bmchk(), which is likely to be 
incorporated within optimization codes.
}
\usage{
   bmchk(par, lower=NULL, upper=NULL, bdmsk=NULL, trace=0, offset=100, shift2bound=TRUE)
}
\arguments{
 \item{par}{A numeric vector of starting values of the optimization function parameters.}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization. Partly for historical reasons, we use the 
     same array during the progress of optimization as an indicator that a 
     parameter is at a lower bound (bdmsk element set to -3) or upper bound (-1).}
 \item{trace}{An integer that controls whether diagnostic information is displayed.
              A positive value displays information, 0 (default) does not.}
 \item{offset}{If provided, is used to detect equality of numbers. That is, two 
     values \code{a} and \code{b} are taken as equal if \code{a + offset} is 
     equal to \code{b + offset}. Default value is 100. Note that in previous versions
     of this code a tolerance \code{tol} was used.}
 \item{shift2bound}{If TRUE, non-masked paramters outside bounds are adjusted 
     to the nearest bound. We then set parchanged = TRUE which implies the 
     original parameters were infeasible.
}
}
\details{
     The bmchk function will check that the bounds exist and are admissible, 
     that is, that there are no lower bounds that exceed upper bounds. 

     There is a check if lower and upper bounds are very close together, in 
     which case a mask is imposed and maskadded is set TRUE. NOTE: it is 
     generally a VERY BAD IDEA to have bounds close together in optimization,
     but here we use a tolerance based on the double precision machine 
     epsilon. Thus it is not a good idea to rely on bmchk() to test if 
     bounds constraints are well-posed.
 
}

\value{
  A list with components:
  \item{bvec}{The vector of parameters, possibly adjusted to bounds. Parameters 
         outside bounds are adjusted to the nearest bound. We let \code{n} be the
         length of this vector.}
  \item{bdmsk}{adjusted input masks}
  \item{bchar}{a set of indicators that give information about the parameters, that is,
               if they are out-of-bounds-low ("-"), at a lower bound ("L"), 
               free ("F"), at an upper bound ("U"), out-of-bounds-high ("+"), 
               masked (fixed) ("M"), of unknown characteristics ("?""), or 
               inadmissible ("!"").}
  \item{lower}{(adjusted) lower bounds. If there are no lower bounds, a vector
         of values equal to -Inf is substituted. Similarly, a single value is
         expanded to a complete vector.
         If any upper and lower bounds are equal (as mediated by \code{offset}),
         we create a mask.}
  \item{upper}{(adjusted) upper bounds}
  \item{nolower}{TRUE if no lower bounds, FALSE otherwise}
  \item{noupper}{TRUE if no upper bounds, FALSE otherwise}
  \item{bounds}{TRUE if there are any bounds, FALSE otherwise}
  \item{admissible}{TRUE if bounds are admissible, FALSE otherwise
         This means no lower bound exceeds an upper bound. That is the bounds 
         themselves are sensible. This condition has nothing to do with the 
         starting parameters.}
  \item{maskadded}{TRUE when a mask has been added because bounds are very close
        or equal, FALSE otherwise. See the code for the implementation.}
  \item{parchanged}{TRUE if parameters are changed by bounds, FALSE otherswise.
        Note that parchanged = TRUE implies the input parameter values were infeasible, 
        that is, violated the bounds constraints.}
  \item{feasible}{TRUE if parameters are within or on bounds, FALSE otherswise.}
  \item{onbound}{TRUE if any parameter is on a bound, FALSE otherswise.
        Note that parchanged = TRUE implies onbound = TRUE, but this is not used inside
        the function. This output value may be important, for example, in using the
        optimization function \code{nmkb} from package \code{dfoptim}.}
}
\examples{
#####################

## cat("25-dimensional box constrained function\n")
## flb <- function(x)
##     { p <- length(x); sum(c(1, rep(4, p-1)) * (x - c(1, x[-p])^2)^2) }

start<-rep(2, 25)
cat("\n start:")
print(start)
lo<-rep(2,25)
cat("\n lo:")
print(lo)
hi<-rep(4,25)
cat("\n hi:")
print(hi)
bt<-bmchk(start, lower=lo, upper=hi, trace=1)
print(bt)

}
\keyword{nonlinear}
\keyword{optimize}
\keyword{upper}
\keyword{lower}
\keyword{bound}
\keyword{mask}
