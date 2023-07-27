multistart <- function(parmat, fn, gr=NULL, lower=-Inf, upper=Inf, 
            method=NULL, hessian=FALSE, control=list(), ...) {
##
   nset <- nrow(parmat)
   npar <- ncol(parmat)
   if (nset < 1) stop("multistart: no starting parameters!")
   ans.ret <- matrix(NA, nrow=nset, ncol=npar+4)
   ans.ret <- data.frame(ans.ret)
   pstring<-colnames(parmat)
   if (is.null(pstring)) {
     pstring <- NULL
     for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
   }  
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence")
   colnames(ans.ret)<-cnames
   row.names(ans.ret)<-1:nset

   for (imeth in 1:nset){
       start <- parmat[imeth, ]
       ans <- optimr(par=start, fn=fn, gr=gr, lower=lower, upper=upper, 
            method=method, hessian=hessian, control=control, ...)
       addvec <- c(ans$par, ans$value, ans$counts[1], ans$counts[2], ans$convergence)
       # print(addvec)
       ans.ret[imeth,] <- addvec
   }
   ans.ret
} ## end of multistart
