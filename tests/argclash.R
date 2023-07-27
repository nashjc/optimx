# argclash.R -- Try to work around dotargs that have a name that clashes 
#     with variable names in function calls
# J C Nash 2021-12-16
# Trying to fix issue raised in 
# https://stackoverflow.com/questions/69033754/maximum-likelihood-estimation-of-three-parameter-reverse-weibull-model-implement/70382556#70382556
# WARNING: It is NOT certain that all name clashes between dotargs and other names
# in functions in optimx will be avoided by the current mechanisms.
# This script does, however, show some tests
rm(list=ls()) # In case we want to ensure clear workspace, delete first #
sqmod<-function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   f<-sum((yy-z)^2)
#   cat("Fv=",f," at ")
#   print(z)
   f
}
sqmod.g <- function(z, x){
   nn<-length(z)
   yy<-x^(1:nn)
   gg<- 2*(z - yy)
}

require(optimx)
require(numDeriv)
sessionInfo()
nn<-2
st<-rep(0.5, nn)
# See if optimx() can handle the clash. 

t2 <- optimx(st, fn=sqmod, x=2)
t2
o1 <- optim(st, fn=sqmod,  x=2)
o1

# Try grchk -- The argclash is fixed in the optimx.run code and in optimr().
tgr <- try(grchk(xpar=st, ffn=sqmod, ggr=sqmod.g, trace=2, x=2))
tgr
# One way that x gets into the dot arguments
xval <- 2
sqmod1 <- function(z){ sqmod(z, x=xval) }
# Another way
dots <- list(x=2)
str(dots)
sqmod2 <- function(z){ sqmod(z, unlist(dots)) }
# simple test
x <- c(1,3)
# Note that x is a listed argument (the 2nd) of numDeriv::grad()
tryg1 <- grad(sqmod1, x )
tryg1
eg1<-sqmod.g(z=x, x=xval)
eg1

cat("sqmod2 for x=2:", sqmod2(x), "\n")
tryg2 <- grad(sqmod2, x )
tryg2

x<-2.0
t2x <- optimx(st, fn=sqmod, x=x)
t2x

t2fm <- optimx(st, fn=sqmod, method="BFGS", x=2)
t2fm

# Following illustrates collision of arguments in gradient after solve
t2fgm <- try(optimx(st, fn=sqmod, gr=sqmod.g, method="BFGS", x=2))
t2fgm


r2 <- optimr(st, fn=sqmod, gr=sqmod.g, method="Rvmmin", control=list(trace=0), x=2)
proptimr(r2)


x<-2.0
r2xn <- optimr(st, fn=sqmod, gr="grfwd", method="Rvmmin", control=list(trace=0), x=x)
proptimr(r2xn)

t2g <- optimx(st, fn=sqmod, gr=sqmod.g, control=list(trace=0), x=2)
t2g
x<-2.0
t2gn <- optimx(st, fn=sqmod,control=list(trace=0), x=x)
t2gn

## Explicit try
xtry<-grchk(xpar=st, ffn=sqmod, ggr=sqmod.g, trace=0, testtol=(.Machine$double.eps^(1/3)), x=x)
xtry

r2g <- optimr(st, fn=sqmod, gr=sqmod.g, method="Rvmmin", control=list(trace=0), x=2)
proptimr(r2g)
x<-2.0
r2xf <- optimr(st, fn=sqmod, gr="grfwd", method="Rvmmin", control=list(trace=0), x=x)
proptimr(r2xf)

