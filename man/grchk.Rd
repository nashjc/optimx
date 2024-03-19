\name{grchk}
\alias{grchk}

\title{Run tests, where possible, on user objective function and (optionally) gradient and hessian}

\description{
   \code{grchk} checks a user-provided R function, \code{ffn}. 
}


\usage{
   grchk(xpar, ffn, ggr, trace=0, testtol=(.Machine$double.eps)^(1/3), ...) 
}


\arguments{
     
    \item{xpar}{
        parameters to the user objective and gradient functions ffn and ggr
        }
    \item{ffn}{
        User-supplied objective function
        }
    \item{ggr}{
        User-supplied gradient function
        }
    \item{trace}{
        set >0 to provide output from grchk to the console, 0 otherwise
        }
    \item{testtol}{
        tolerance for equality tests
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}


\details{
\tabular{ll}{
Package: \tab grchk\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
\code{numDeriv} is used to numerically approximate the gradient of function \code{ffn}
and compare this to the result of function \code{ggr}.

}


\value{
\code{grchk} returns a single object \code{gradOK} which is TRUE if the differences 
between analytic and approximated gradient are small as measured by the tolerance 
\code{testtol}.

This has attributes "ga" and "gn" for the analytic and numerically approximated gradients,
and "maxdiff" for the maximum absolute difference between these vectors.

At the time of preparation, there are no checks for validity of the gradient code in
\code{ggr} as in the function \code{fnchk}.


}


\author{

    John C. Nash
}

\examples{
# Would like examples of success and failure. What about "near misses"?
cat("Show how grchk works\n")
require(numDeriv)
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

jonesg2 <- function(xx) {
  gx <- 1
  gy <- 2
  gg <- - c(gx, gy)
}


xx <- c(1, 2)

gcans <- grchk(xx, jones, jonesg, trace=1, testtol=(.Machine$double.eps)^(1/3))
gcans

gcans2 <- grchk(xx, jones, jonesg2, trace=1, testtol=(.Machine$double.eps)^(1/3))
gcans2




}

\keyword{optimize}

