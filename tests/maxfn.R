# maxfn.R
##  author: John C. Nash

# solution: 1:n has maximum of 10.

maxfn<-function(x) {# fn to be MAXIMIZED
  # max = 10 at 1:n
  n<-length(x)
  ss<-seq(1,n)
  f<-10-sum((x-ss)^2)
  f
}

maxfn.g <- function(x) { # gradient
   n <- length(x)
   ss<-seq(1,n)
   gg <- -2*(x-ss)
   gg
}

maxfn.h <- function(x) { # gradient
  n <- length(x)
  hh<-rep(-2, n)
  hh <- diag(hh)
  hh
}

# solution: 1:n has minimum of -10.

negmaxfn<-function(x) {# explicit negative of maxfn
  f<-(-1)*maxfn(x)
  return(f)
}
negmaxfn.g<-function(x) {# explicit negative of maxfn
  gg<-(-1)*maxfn.g(x)
  gg
}
negmaxfn.h<-function(x) {# explicit negative of maxfn
  hh<-(-1)*maxfn.h(x)
  hh
}

