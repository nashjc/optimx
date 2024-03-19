# Simple Test Function 1:
simfun.f = function(x) {
  fun <- sum(x^2 )
  ## if (trace) ... to be fixed
  print(c(x = x, fun = fun))
  fun
}
simfun.g = function(x) {
  grad<-2.0*x
  grad
}
simfun.h = function(x) {
  n<-length(x)
  t<-rep(2.0,n)
  hess<-diag(t)
}

fghfn <- function(x){ # dotargs removed 230619
  f <- simfun.f(x) # dotargs removed 230619
  g <- simfun.g(x)
  h <- simfun.h(x)
                          
  attr(f,"gradient") <- g
  attr(f,"hessian") <- h
  f
}

library(optimx)

strt <- c(1,2,3)
ansfgh <- optimr(strt, simfun.f, simfun.g, simfun.h, method="nlm",
                 hessian=TRUE, control=list(trace=2))
cat("ansfgh$counts:"); print(ansfgh$counts)

proptimr(ansfgh) # compact output of result
cat("\n All result object:")
print(ansfgh)
anlm<-nlm(fghfn,strt, hessian=TRUE)
print(anlm)          