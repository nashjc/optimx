\name{grback}
\alias{grback}

\title{Backward difference numerical gradient approximation.}

\description{
   \code{grback} computes the backward difference approximation to the gradient of 
          user function \code{userfn}. 
}


\usage{
   grback(par, userfn, fbase=NULL, env=optsp, ...)
}


\arguments{
     
    \item{par}{
        parameters to the user objective function userfn
        }
    \item{userfn}{
        User-supplied objective function
        }
    \item{fbase}{
        The value of the function at the parameters, else NULL. This is to save
        recomputing the function at this point.
        }
    \item{env}{
        Environment for scratchpad items (like \code{deps} for approximation 
        control in this routine). Default \code{optsp}.
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}


\details{
\tabular{ll}{
Package: \tab grback\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
}


\value{
\code{grback} returns a single vector object \code{df} which approximates the 
gradient of userfn at the parameters par. The approximation is controlled by a
global value \code{optderiveps} that is set when the package is attached.

}


\author{

    John C. Nash
}

\examples{
cat("Example of use of grback\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}

xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grback(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to analytic gradient:\n")
print(ga)

cat("change the step parameter to 1e-4\n")
optsp$deps <- 1e-4
gn2<-grback(xx,myfn, shift=0)
print(gn2)

}

\keyword{optimize}
