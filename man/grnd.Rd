\name{grnd}
\alias{grnd}
\encoding{UTF-8}
\title{A reorganization of the call to numDeriv grad() function.}
\concept{minimization}
\concept{maximization}
\description{
        Provides a wrapper for the numDeriv approximation to the
        gradient of a user supplied objective function \code{userfn}.
}
\usage{
       grnd(par, userfn, ...)
}
\arguments{
 \item{par}{A vector of parameters to the user-supplied function \code{fn}}
 \item{userfn}{A user-supplied function }
 \item{...}{Other data needed to evaluate the user function.}
}
\details{
  The Richardson method is used in this routine.
}
\value{
  \code{grnd} returns an approximation to the gradient of the function userfn
}
\examples{
cat("Example of use of grnd\n")
require(numDeriv)
myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grnd(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)
}
\keyword{nonlinear}
\keyword{optimize}

