\name{polyopt}
\alias{polyopt}
\encoding{UTF-8}
\title{General-purpose optimization - sequential application of methods}
\concept{minimization}
\concept{maximization}
\description{
  Multiple minimization methods are applied in sequence to a single problem,
  with the output parameters of one method being used to start the next.
}
\usage{
polyopt(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            methcontrol=NULL, hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{methcontrol}{An data frame of which each row gives an optimization method, a maximum
       number of iterations and a maximum number of function evaluations allowed for that
       method. Each method will be executed in turn until either the maximum iterations
       or function evaluations are completed, whichever is first. The next method is then
       executed starting with the best parameters found so far, else the function exits.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim().}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  See \code{optimr()} for other details.
  
  Note that this function does not (yet?) make use of a hess function
  to compute the hessian.
}
\value{

An array with one row per method. Each row contains:

   \item{par}{The best set of parameters found for the method in question.}
   \item{value}{The value of ‘fn’ corresponding to ‘par’.}
   \item{counts}{ A two-element integer vector giving the number of calls to
          ‘fn’ and ‘gr’ respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to ‘fn’
          to compute a finite-difference approximation to the gradient.}
   \item{convergence}{ An integer code. ‘0’ indicates successful completion}
   \item{ message}{ A character string giving any additional information returned
          by the optimizer, or ‘NULL’.}
   \item{hessian}{ Always NULL for this routine.}
}
\source{
See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.
}
\examples{
fnR <- function (x, gs=100.0) 
{
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n - 1)]
    sum(gs * (x1 - x2^2)^2 + (1 - x2)^2)
}
grR <- function (x, gs=100.0) 
{
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2 * (x[1] - 1) + 4*gs * x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n - 1)
        g[ii] <- 2 * (x[ii] - 1) + 4 * gs * x[ii] * (x[ii]^2 - x[ii + 
            1]) + 2 * gs * (x[ii] - x[ii - 1]^2)
    }
    g[n] <- 2 * gs * (x[n] - x[n - 1]^2)
    g
}

x0 <- rep(pi, 4)
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(1000, 100), maxfeval= c(1000, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans
mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(100, 100), maxfeval= c(100, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans

mc <- data.frame(method=c("Nelder-Mead","Rvmmin"), maxit=c(10, 100), maxfeval= c(10, 1000))

ans <- polyopt(x0, fnR, grR, methcontrol=mc, control=list(trace=0))
ans



}
\keyword{nonlinear}
\keyword{optimize}
