\name{fnchk}
\alias{fnchk}
\title{Run tests, where possible, on user objective function}
\description{
   \code{fnchk} checks a user-provided R function, \code{ffn}. 
}
\usage{
   fnchk(xpar, ffn, trace=0, ... )
}
\arguments{
     
    \item{xpar}{
        the (double) vector of parameters to the objective funcion
        }
    \item{ffn}{
        a user-provided function to compute the objective function
        }
    \item{trace}{
	set >0 to provide output from fnchk to the console, 0 otherwise
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
}
\details{
\code{fnchk} attempts to discover various errors in function setup in user-supplied
functions primarily intended for use in optimization calculations. There are always 
more conditions that could be tested!    
}
\value{
  The output is a list consisting of 
list(fval=fval, infeasible=infeasible, excode=excode, msg=msg)
  \item{fval}{The calculated value of the function at parameters \code{xpar} if the function
      can be evaluated.}
  \item{infeasible}{FALSE if the function can be evaluated, TRUE if not.}
  \item{excode}{An exit code, which has a relationship to }
  \item{msg}{A text string giving information about the result of the function check: Messages and
    the corresponding values of \code{excode} are:
    \itemize{
        \item{fnchk OK;}{ \code{excode} = 0; 
                       \code{infeasible} = FALSE}
        \item{Function returns INADMISSIBLE;}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
        \item{Function returns a vector not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns a list not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns a matrix list not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returns an array not a scalar;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returned not length 1, despite not vector, matrix or array;}
                    { \code{excode} = -4; \code{infeasible} = TRUE}
        \item{Function returned non-numeric value; \code{excode} = 0;}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
        \item{Function returned Inf or NA (non-computable);}
                    { \code{excode} = -1; \code{infeasible} = TRUE}
    }
  }
}


\author{
    John C. Nash <nashjc@uottawa.ca>
}

\examples{
# Want to illustrate each case.
# Ben Bolker idea for a function that is NOT scalar
# rm(list=ls())
# library(optimx)
sessionInfo()
benbad<-function(x, y){
  # y may be provided with different structures
  f<-(x-y)^2
} # very simple, but ...

y<-1:10
x<-c(1)
cat("fc01: test benbad() with y=1:10, x=c(1)\n")
fc01<-fnchk(x, benbad, trace=4, y)
print(fc01)

y<-as.vector(y)
cat("fc02: test benbad() with y=as.vector(1:10), x=c(1)\n")
fc02<-fnchk(x, benbad, trace=1, y)
print(fc02)

y<-as.matrix(y)
cat("fc03: test benbad() with y=as.matrix(1:10), x=c(1)\n")
fc03<-fnchk(x, benbad, trace=1, y)
print(fc03)

y<-as.array(y)
cat("fc04: test benbad() with y=as.array(1:10), x=c(1)\n")
fc04<-fnchk(x, benbad, trace=1, y)
print(fc04)

y<-"This is a string"
cat("test benbad() with y a string, x=c(1)\n")
fc05<-fnchk(x, benbad, trace=1, y)
print(fc05)

cat("fnchk with Rosenbrock\n")
fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
xtrad<-c(-1.2,1)
ros1<-fnchk(xtrad, fr, trace=1)
print(ros1)
npar<-2
opros<-list2env(list(fn=fr, gr=NULL, hess=NULL, MAXIMIZE=FALSE, PARSCALE=rep(1,npar), FNSCALE=1,
                     KFN=0, KGR=0, KHESS=0, dots=NULL))
uros1<-fnchk(xtrad, fr, trace=1)
print(uros1)


}
\keyword{optimize}

