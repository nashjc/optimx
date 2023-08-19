% File src/library/stats/man/getresult.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{gHgenb}
\alias{gHgenb}
\encoding{UTF-8}
\title{Generate gradient and Hessian for a function at given parameters.}
\description{
  \code{gHgenb} is used to generate the gradient and Hessian of an objective
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
  gHgenb(par, fn, gr=NULL, hess=NULL, bdmsk=NULL, lower=NULL, upper=NULL,
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
 \item{bdmsk}{An integer vector of the same length as \code{par}. When an element
     of this vector is 0, the corresponding parameter value is fixed (masked) 
     during an optimization. Non-zero values indicate a parameter is free (1),
     at a lower bound (-3) or at an upper bound (-1), but this routine only
     uses 0 values.}
 \item{lower}{Lower bounds for parameters in \code{par}.}
 \item{upper}{Upper bounds for parameters in \code{par}.}
 \item{control}{A list of controls to the function. Currently 
                  asymptol (default of 1.0e-7 which tests for asymmetry of Hessian approximation
                  (see code for details of the test); 
                  ktrace, a logical flag which, if TRUE, monitors the progress 
                  of gHgenb (default FALSE), and 
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
      \item{\code{nbm}}{  The number of active bounds and masks.}
  }
}
\examples{
require(numDeriv)
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
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}


maxfn<-function(x, top=10) {
      	n<-length(x)
	ss<-seq(1,n)
	f<-top-(crossprod(x-ss))^2
	f<-as.numeric(f)
	return(f)
}

negmaxfn<-function(x) {
	f<-(-1)*maxfn(x)
	return(f)
}

parx<-rep(1,4)
lower<-rep(-10,4)
upper<-rep(10,4)
bdmsk<-c(1,1,0,1) # masked parameter 3
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgenb(parx,genrose.f)
cat("results of gHgenb for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgenb(parx,genrose.f,genrose.g)
cat("results of gHgenb for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("*****************************************\n")
parx<-rep(0.9,4)
fval<-genrose.f(parx)
gval<-genrose.g(parx)
Ahess<-genrose.h(parx)
gennog<-gHgenb(parx,genrose.f,control=list(ktrace=TRUE), gs=9.4)
cat("results of gHgenb with gs=",9.4," for genrose without gradient code at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
cat("\n\n")
geng<-gHgenb(parx,genrose.f,genrose.g, control=list(ktrace=TRUE))
cat("results of gHgenb for genrose at ")
print(parx)
print(gennog)
cat("compare to g =")
print(gval)
cat("and Hess\n")
print(Ahess)
gst<-5
cat("\n\nTest with full calling sequence and gs=",gst,"\n")
gengall<-gHgenb(parx,genrose.f,genrose.g,genrose.h, control=list(ktrace=TRUE),gs=gst)
print(gengall)


top<-25
x0<-rep(2,4)
cat("\n\nTest for maximization and top=",top,"\n")
cat("Gradient and Hessian will have sign inverted")
maxt<-gHgen(x0, maxfn, control=list(ktrace=TRUE), top=top)
print(maxt)

cat("test against negmaxfn\n")
gneg <- grad(negmaxfn, x0)
Hneg<-hessian(negmaxfn, x0)
# gdiff<-max(abs(gneg-maxt$gn))/max(abs(maxt$gn))
# Hdiff<-max(abs(Hneg-maxt$Hn))/max(abs(maxt$Hn))
# explicitly change sign 
gdiff<-max(abs(gneg-(-1)*maxt$gn))/max(abs(maxt$gn))
Hdiff<-max(abs(Hneg-(-1)*maxt$Hn))/max(abs(maxt$Hn))
cat("gdiff = ",gdiff,"  Hdiff=",Hdiff,"\n")



}
\keyword{nonlinear}
\keyword{optimize}
