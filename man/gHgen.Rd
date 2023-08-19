% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{gHgen}
\alias{gHgen}
\encoding{UTF-8}
\title{Generate gradient and Hessian for a function at given parameters.}
\description{
  \code{gHgen} is used to generate the gradient and Hessian of an objective
  function used for optimization. If a user-provided gradient function 
  \code{gr} is available it is used to compute the gradient, otherwise 
  package \code{numDeriv} is used. If a user-provided Hessian function
  \code{hess} is available, it is used to compute a Hessian. Otherwise, if
  \code{gr} is available, we use the function \code{jacobian()} from
  package \code{numDeriv} to compute the Hessian. In both these cases we
  check for symmetry of the Hessian. Computational Hessians are commonly
  NOT symmetric. If only the objective function \code{fn} is provided, then
  the Hessian is approximated with the function \code{hessian} from 
  package \code{numDeriv} which guarantees a symmetric matrix. 
}
\usage{
  gHgen(par, fn, gr=NULL, hess=NULL,
      control=list(ktrace=0), ...) 
}
\arguments{
 \item{par}{Set of parameters, assumed to be at a minimum of the function \code{fn}.}
 \item{fn}{Name of the objective function.}
 \item{gr}{(Optional) function to compute the gradient of the objective function. If present,
     we use the Jacobian of the gradient as the Hessian and avoid one layer of numerical
     approximation to the Hessian.} 
 \item{hess}{(Optional) function to compute the Hessian of the objective function. This
     is rarely available, but is included for completeness.}
 \item{control}{A list of controls to the function. Currently 
                  asymptol (default of 1.0e-7 which tests for asymmetry of Hessian approximation
                  (see code for details of the test); 
                  ktrace, a logical flag which, if TRUE, monitors the progress 
                  of gHgen (default FALSE), and 
                  stoponerror, defaulting to FALSE to NOT stop when there is
                  an error or asymmetry of Hessian. Set TRUE to stop.}
 \item{\dots}{Extra data needed to compute the function, gradient and Hessian.}
}

\details{
  None
}
\value{
  \code{ansout} a list of four items, 
  \itemize{
      \item{\code{gn}}{  The approximation to the gradient vector.}
      \item{\code{Hn}}{  The approximation to the Hessian matrix.}
      \item{\code{gradOK}}{  TRUE if the gradient has been computed acceptably. FALSE otherwise.}
      \item{\code{hessOK}}{  TRUE if the gradient has been computed acceptably and passes the
             symmetry test. FALSE otherwise.}
      \item{\code{nbm}}{  Always 0. The number of active bounds and masks.
             Present to make function consistent with \code{gHgenb}.}
  }
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
ans100fgh<-  gHgen(trad, genrose.f, gr=genrose.g, hess=genrose.h,
      control=list(ktrace=1)) 
print(ans100fgh)
ans100fg<-  gHgen(trad, genrose.f, gr=genrose.g, 
      control=list(ktrace=1)) 
print(ans100fg)
ans100f<-  gHgen(trad, genrose.f, control=list(ktrace=1)) 
print(ans100f)
ans10fgh<-   gHgen(trad, genrose.f, gr=genrose.g, hess=genrose.h,
      control=list(ktrace=1), gs=10) 
print(ans10fgh)
ans10fg<-   gHgen(trad, genrose.f, gr=genrose.g, 
      control=list(ktrace=1), gs=10) 
print(ans10fg)
ans10f<-   gHgen(trad, genrose.f, control=list(ktrace=1), gs=10) 
print(ans10f)

}

\keyword{nonlinear}
\keyword{optimize}
