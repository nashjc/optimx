\name{optimr}
\alias{optimr}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optimr} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end, in fact using the
  calling sequence of the R function \code{optim()}.
}
\usage{
optimr(par, fn, gr=NULL, hess=NULL, method=NULL, lower=-Inf, upper=Inf, 
          hessian=FALSE, control=list(), ...)
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

   If \code{gr} is \code{NULL}, then this is passed forward and 
   whatever default action is specified for the 
   chosen method for the case of a null gradient code
   is used. For many methods, this is a finite-difference approximation, 
   but some methods require user input for the gradient and will fail
   otherwise. In such cases, we try to return \code{convergence} of 9998. 
   
   If \code{gr} is a character string, then that string is taken as the name of 
   a gradient approximation function, for example, "grfwd", "grback" and
   "grcentral" for standard forward, backward and central approximations.
   Method "grnd" uses the \code{grad()} function from package \code{numDeriv}.
 }   
  \item{hess}{A function to return (as a matrix) the hessian for those methods that
   can use this information.

  If \code{hess} is the character string "approx", then ??

  If \code{hess} is \code{NULL}, then this is passed forward and 
   whatever default action is specified for the 
   chosen method for the case of a null gradient code is used. 
 }
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints. A small set of methods can handle masks, that is,
   fixed parameters, and these can be specified by making the lower and upper bounds
   equal to the starting value. (It is possible that the starting value could be different
   from the lower/upper bounds set,
   but this behaviour has NOT yet been defined and users are cautioned.)}
 \item{method}{A character string giving the name of the optimization method to be
    applied. See the list \code{allmeth} in file 
    \code{ctrldefault.R} which is part of this package.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. Note that this will NOT necessarily
       use the same approximation as may be provided by the method called. Instead, 
       the function \code{hessian()} from package \code{numDeriv} is used if no gradient
       \code{gr} is supplied, else the function \code{jacobian()} from \code{numDeriv}
       is applied to the gradient function \code{gr}.}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{Further arguments to be passed to \code{fn} 
    and \code{gr} if needed for computation of these quantities; otherwise, further 
    arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} should be matched exactly.

  By default \code{optimr} performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this, and this
  control can be used with \code{optimr} but is deprecated. Moreover, 
  if \code{control$maximize} is set, it will
  take precedence over \code{control$fnscale}. Generally it is a BAD IDEA
  to use both mechanisms simultaneously. 

  Possible method choices are specified by the list \code{allmeth} in the file
  \code{ctrldefault.R} which is part of this package. 

  If no method is specified, the method specified by \code{defmethod} in file 
  \code{ctrldefault.R} (which is part of this package) will be attempted. 
  
  Function \code{fn} must return a finite scalar value at the initial set
  of parameters. Some methods can handle a returned value \code{NA} or \code{Inf} 
  if the function cannot be evaluated at the supplied value. However, other 
  methods, of which \code{"L-BFGS-B"} is known to be a case, require that 
  the values returned should always be finite. It is recommended that user functions
  ALWAYS return a usable value. Note that the control \code{badval} in 
  \code{ctrldefault.R} give a possible number that could be returned.

  For details of methods, please consult the documentation of the individual methods.
  (The NAMESPACE file lists the packages from which functions are imported.)
  Note that method \code{"hjn"} is a conservative implementation of a Hooke and 
  Jeeves (1961) and is part of this package. It is provided as a simple example of
  a very crude optimization method; it is NOT intended as a production method, but
  may be useful for didactic purposes.

  The \code{control} argument is a list that can supply any of the
  components in the file \code{ctrldefault.R} which is part of this 
  package. It may supply controls that are
  useful or required for particular methods, but users are warned to be careful to
  ensure that extraneous or incorrect components and values are not passed.
  Some \code{control} elements apply only to some methods. 
  See individual packages for details. \code{optimr} does not support all the
  possible controls for all methods.
  
  A particular case is the method "bobyqa", where the control \code{rhobeg=0} 
  gives a set of controls that depend on the bounds supplied. This choice is
  only in the current package. Unspecified or negative control \code{rhobeg=0}
  gives the minqa defaults. Positive value of this control (and optionally 
  control \code{rhoend}) supply those values. 
  See inst/doc/examples/specctrlhobbs.R.
  
  
  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Apparently no other attributes of \code{par}
  are copied over, but this may need to be verified, especially if parameters
  are passed to non-R routines.
  
  CAUTION: because there is a seldom-used parameter \code{hess}, you should NOT
  make a call like 
  
ans <- optimr(start, myf, myg, lower, upper)
     
  or you will likely get wrong results. Instead use
  
ans <-  optimr(start, myf, myg, lower=lower, upper=upper)
}

\value{

    A list with components:
  \describe{
   \item{par}{The best set of parameters found.}
   \item{value}{The value of ‘fn’ corresponding to ‘par’.}
   \item{counts}{ A two-element integer vector giving the number of calls to
          ‘fn’ and ‘gr’ respectively. This excludes those calls needed
          to compute the Hessian, if requested, and any calls to ‘fn’
          to compute a finite-difference approximation to the gradient.
          NOT available to be reported for some methods, e.g.,  \code{lbfgs}.}
   \item{convergence}{ An integer code. ‘0’ indicates successful completion. The
   documentation for function \code{opm()} gives some other possible values and
   their meaning. }
   \item{ message}{ A character string giving any additional information returned
          by the optimizer, or ‘NULL’.}
   \item{hessian}{If requested, an approximation to the hessian of ‘fn’
        at the final parameters.}
  }
}
\examples{
 # Simple Test Function 1:
simfun.f = function(x) {
     fun <- sum(x^2 )
## if (trace) ... to be fixed
	print(c(x = x, fun = fun))
     fun
}
simfun.g = function(x) {
     grad<-2.0*x
     grad
}
simfun.h = function(x) {
     n<-length(x)
     t<-rep(2.0,n)
     hess<-diag(t)
}

strt <- c(1,2,3)
ansfgh <- optimr(strt, simfun.f, simfun.g, simfun.h, method="nlm",
     hessian=TRUE, control=list(trace=2))
proptimr(ansfgh) # compact output of result


}

\references{

See the manual pages for \code{optim()}.

 Hooke R. and Jeeves, TA (1961). Direct search solution of numerical and statistical problems. 
   Journal of the Association for Computing Machinery (ACM). 8 (2): 212–229.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

 Nocedal J, and Wright SJ (1999). Numerical optimization. New York: Springer. 2nd Edition 2006.
 
}
\keyword{nonlinear}
\keyword{optimize}
