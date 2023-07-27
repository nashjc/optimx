% File src/library/stats/man/kktchk.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later


\name{kktchk}
\alias{kktchk}
\encoding{UTF-8}
\title{Check Kuhn Karush Tucker conditions for a supposed function minimum}
\concept{minimization}
\concept{maximization}
\description{
  Provide a check on Kuhn-Karush-Tucker conditions based on quantities
  already computed. Some of these used only for reporting.
}
\usage{
     kktchk(par, fn, gr, hess=NULL, upper=NULL, lower=NULL, 
                 maximize=FALSE, control=list(), ...)
}
\arguments{
 \item{par}{A vector of values for the parameters which are supposedly optimal.}
 \item{fn}{The objective function}
 \item{gr}{The gradient function}
 \item{hess}{The Hessian function}
 \item{upper}{Upper bounds on the parameters}
 \item{lower}{Lower bounds on the parameters}
 \item{maximize}{Logical TRUE if function is being maximized. Default FALSE.}
 \item{control}{A list of controls for the function}
 \item{...}{The dot arguments needed for evaluating the function and gradient and hessian}

}
\details{
  kktchk computes the gradient and Hessian measures for BOTH unconstrained and 
  bounds (and masks) constrained parameters, but the kkt measures are evaluated
  only for the constrained case.

  Note that evaluated Hessians are often not symmetric, and many, possibly most,
  examples will fail the \code{is.Symmetric()} function. In such cases, the 
  check on the Hessian uses the mean of the Hessian and its transpose.

}

\value{
  The output is a list consisting of 
  \item{gmax}{The absolute value of the largest gradient component in magnitude.}
  \item{evratio}{The ratio of the smallest to largest Hessian eigenvalue. Note that this
     may be negative.}
  \item{kkt1}{A logical value that is TRUE if we consider the first (i.e., gradient) 
     KKT condition to be satisfied. WARNING: The decision is dependent on tolerances and
     scaling that may be inappropriate for some problems.}
  \item{kkt2}{A logical value that is TRUE if we consider the second (i.e., positive
      definite Hessian) KKT condition to be satisfied. WARNING: The decision is dependent 
      on tolerances and scaling that may be inappropriate for some problems.}
  \item{hev}{The calculated hessian eigenvalues, sorted largest to smallest. 
     Sorting is a property of the \code{eigen()} function.}
  \item{ngatend}{The computed (unconstrained) gradient at the solution parameters.}
  \item{nnatend}{The computed (unconstrained) hessian at the solution parameters.}
}
\seealso{

    \code{\link{optim}}

}
\examples{
cat("Show how kktc works\n")

# require(optimx)

jones<-function(xx){
  x<-xx[1]
  y<-xx[2]
  ff<-sin(x*x/2 - y*y/4)*cos(2*x-exp(y))
  ff<- -ff
}

jonesg <- function(xx) {
  x<-xx[1]
  y<-xx[2]
  gx <-  cos(x * x/2 - y * y/4) * ((x + x)/2) * cos(2 * x - exp(y)) - 
    sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * 2)
  gy <- sin(x * x/2 - y * y/4) * (sin(2 * x - exp(y)) * exp(y)) - cos(x * 
             x/2 - y * y/4) * ((y + y)/4) * cos(2 * x - exp(y))
  gg <- - c(gx, gy)
}

ans <- list() # to ensure structure available
# If optimx package available, the following can be run.
# xx<-0.5*c(pi,pi)
# ans <- optimr(xx, jones, jonesg, method="Rvmmin")
# ans

ans$par <- c(3.154083, -3.689620)

kkans <- kktchk(ans$par, jones, jonesg)
kkans



}
\keyword{nonlinear}
\keyword{optimize}
