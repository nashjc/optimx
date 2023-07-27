\name{hesschk}
\alias{hesschk}
\title{Run tests, where possible, on user objective function and (optionally) gradient and hessian}
\description{
   \code{hesschk} checks a user-provided R function, \code{ffn}. 
}
\usage{
   hesschk(xpar, ffn, ggr, hhess, trace=0, testtol=(.Machine$double.eps)^(1/3), ...)
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
    \item{hhess}{
        User-supplied Hessian function
        }
    \item{trace}{
        set >0 to provide output from hesschk to the console, 0 otherwise
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
Package: \tab hesschk\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
\code{numDeriv} is used to compute a numerical approximation to the Hessian
matrix. If there is no analytic gradient, then the \code{hessian()} function 
from \code{numDeriv} is applied to the user function \code{ffn}. Otherwise, 
the \code{jacobian()} function of \code{numDeriv} is applied to the \code{ggr}
function so that only one level of differencing is used.
}
\value{
The function returns a single object \code{hessOK} which is TRUE if the 
analytic Hessian code returns a Hessian matrix that is "close" to the 
numerical approximation obtained via \code{numDeriv}; FALSE otherwise.

\code{hessOK} is returned with the following attributes:
\itemize{
     \item{"nullhess"}{Set TRUE if the user does not supply a function to compute the Hessian.}
     \item{"asym"}{Set TRUE if the Hessian does not satisfy symmetry conditions to
           within a tolerance. See the \code{hesschk} for details.}
     \item{"ha"}{The analytic Hessian computed at paramters \code{xpar} using \code{hhess}.}
     \item{"hn"}{The numerical approximation to the Hessian computed at paramters \code{xpar}.}
     \item{"msg"}{A text comment on the outcome of the tests.}
}
}
\author{
    John C. Nash
}
\examples{
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
#		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

trad<-c(-1.2,1)
ans100<-hesschk(trad, genrose.f, genrose.g, genrose.h, trace=1)
print(ans100)
ans10<-hesschk(trad, genrose.f, genrose.g, genrose.h, trace=1, gs=10)
print(ans10)


}
\keyword{optimize}

