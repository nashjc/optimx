\name{grfwd}
\alias{grfwd}
\alias{optsp}

\title{Forward difference numerical gradient approximation.}

\description{
   \code{grfwd} computes the forward difference approximation to the gradient of 
          user function \code{userfn}. 
}
\usage{
   grfwd(par, userfn, fbase=NULL, env=optsp, ...)
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
Package: \tab grfwd\cr
Depends: \tab R (>= 2.6.1)\cr
License: \tab GPL Version 2.\cr
}  
}

\value{
\code{grfwd} returns a single vector object \code{df} which approximates the 
gradient of userfn at the parameters par. The approximation is controlled by a
global value \code{optderiveps} that is set when the package is attached.
}
\author{
    John C. Nash
}

\examples{
cat("Example of use of grfwd\n")

myfn<-function(xx, shift=100){
    ii<-1:length(xx)
    result<-shift+sum(xx^ii)
}
xx<-c(1,2,3,4)
ii<-1:length(xx)
print(xx)
gn<-grfwd(xx,myfn, shift=0)
print(gn)
ga<-ii*xx^(ii-1)
cat("compare to\n")
print(ga)
}

\keyword{optimize}

