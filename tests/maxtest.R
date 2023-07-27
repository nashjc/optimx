# maxtest.R
##  author: John C. Nash
# rm(list=ls())
cat("maxfn.R -- test that maximize=TRUE works correctly\n")

##?? Notes Dec 2 2021

##?? When we are using grfwd, the sign of the gradient is wrong. Why??


require(optimx)
# source("optimx/tests/maxfn.R")
source("maxfn.R")
sessionInfo()


n<-4
xx<-rep(1,n) # start at all 1s

rv0 <- Rvmmin(xx, maxfn, maxfn.g, control=list(maximize=TRUE, trace=2))
rv0
rv0f <- Rvmmin(xx, maxfn, "grfwd", control=list(maximize=TRUE, trace=2))
rv0f

nrv0 <- Rvmmin(xx, negmaxfn, negmaxfn.g, control=list(trace=4))
nrv0

nrv0f <- Rvmmin(xx, negmaxfn, "grfwd", control=list(trace=4))
nrv0f


# Conflicting controls -- 'maximize' takes precedence over 'fnscale'
ansconf0<-optimr(xx,maxfn, gr=maxfn.g, method="Rvmmin", control=list(maximize=TRUE, trace=2))
proptimr(ansconf0) # should work OK

ansconf0f<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE, trace=4))
proptimr(ansconf0f) # should work OK
ansconf0c<-optimr(xx,maxfn, gr="grcentral", method="Rvmmin", control=list(maximize=TRUE, trace=2))
proptimr(ansconf0c) # should work OK
ansconf0n<-optimr(xx,maxfn, gr="grnd", method="Rvmmin", control=list(maximize=TRUE, trace=2))
proptimr(ansconf0n) # should work OK
ansconf1<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE, fnscale=1))
proptimr(ansconf1) # should work OK
# Warning
ansconf2<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=FALSE, fnscale=-1))
proptimr(ansconf2) # should work OK


maxall <- opm(xx, maxfn, gr=maxfn.g, hess=maxfn.h, method="ALL", control=list(maximize=TRUE))
summary(maxall, order=value)

# specifying both maximize and fnscale gives maximize precedence with no message
ansmaxboth<-optimr(xx,maxfn, gr="grfwd", method="Rvmmin", control=list(maximize=TRUE, fnscale=-1.0, trace=2))
proptimr(ansmaxboth)
# 
# cat("using the negmax function should give same parameters\n")
ansnegmax<-optimr(xx,negmaxfn, gr="grfwd",  method="Rvmmin", control=list(trace=0))
proptimr(ansnegmax)# function value should be -10 however
# 
ansmaxnm<-optimr(xx,maxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=TRUE,trace=1, maxit=2000))
proptimr(ansmaxnm)# try with Nelder-Mead
# 
ansmaxnmbad<-optimr(xx,negmaxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=TRUE,trace=1))
# THIS SHOULD FAIL UNBOUNDED
proptimr(ansmaxnmbad)
# 
ansmaxnmbad2<-optimr(xx,negmaxfn, gr="grfwd", method="Nelder-Mead", control=list(maximize=FALSE,
                             fnscale=-1, trace=2))
# THIS SHOULD FAIL UNBOUNDED ALSO
proptimr(ansmaxnmbad2)
# 
# 
# # Test Carlo Lapid suggested fix for optimr()  180313
# 
# amaxo<-optimr(xx, maxfn, method="L-BFGS-B", control=list(maximize=TRUE, trace=0))
# proptimr(amaxo)
# 
# cat("using the negmax function should give same parameters\n")
# anegmaxo<-optimr(xx,negmaxfn, method="L-BFGS-B", control=list(trace=0))
# proptimr(anegmaxo)
# 
# 
# 
# cat("WARNING -- this example should FAIL\n")
# cat("maximize=TRUE is NOT set up in hjn()\n")
# # 160706 -- not set up to maximize, except through optimr perhaps
# n<-4
# xx<-rep(1,n)
# ansmax<-try(hjn(xx,maxfn, control=list(maximize=TRUE,trace=1, maxit=10, maxfeval=2000)))
# print(ansmax)
# 
cat("\nTry to maximize through optimr()\n")
anshjno<-try(optimr(xx,maxfn, method="hjn", control=list(maximize=TRUE,trace=1, maxit=10, maxfeval=2000)))
proptimr(anshjno)
# Check kkt conditions
xxh<-anshjno$par
xxhkkt<-kktchk(xxh, maxfn, "grcentral", maximize=TRUE, control=list(trace=2))
