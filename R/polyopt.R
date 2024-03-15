polyopt <- function(par, fn, gr=NULL, lower=-Inf, upper=Inf, 
            methcontrol=NULL, hessian=FALSE, control=list(), ...) {
   nmeth <- nrow(methcontrol)
   npar <- length(par)
   if (nmeth < 1) stop("polyopt: no starting parameters!")
   ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+4)
   ans.ret <- data.frame(ans.ret)
   pstring<-names(par)
   if (is.null(pstring)) {
     pstring <- NULL
     for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
   }  
   cnames <- c(pstring, "value", "fevals", "gevals", "convergence")
   colnames(ans.ret) <- cnames
   row.names(ans.ret) <- 1:nmeth
   cpar <- par # use initial parameters here
# We want a lot of checks here
   for (imeth in 1:nmeth){
       meth <- as.character(methcontrol[imeth,1]) # name of method
       cat("Method ",imeth," :",meth,"\n")
       control$maxit <- methcontrol[[imeth,2]]
       control$maxfeval <- methcontrol[imeth,3]
## Do we need this -- should be able to handle infinite or null bounds?
       if ((is.null(lower) && is.null(upper)) || (is.infinite(lower) && is.infinite(upper))) {
#           cat("calling optimr unconstrained\n")
           ans <- optimr(cpar, fn=fn, gr=gr, 
            method=meth, hessian=FALSE, control=control, ...)
       } else {
#          cat("calling optimr bounded\n")
          ans <- optimr(cpar, fn=fn, gr=gr, lower=lower, upper=upper, 
            method=meth, hessian=FALSE, control=control, ...)
       }
#       ans$xtimes <- NA # 160703 -- not yet available
       addvec <- c(ans$par, ans$value, ans$counts[1], ans$counts[2], 
                     ans$convergence)
#       print(addvec)
       cpar <- ans$par # copy the parameters for next method
#       print(cpar)
       ans.ret[imeth,] <- addvec
   }
   ans.ret
} ## end of polyopt
