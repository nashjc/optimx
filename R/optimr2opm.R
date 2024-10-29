# optimr2opm function to put optimr solutions into opm form
# par, value, counts, convergence
# do we want full optimr form (2 element counts, message and hessian), even
#  if we make those emtpy but they match structure?
#
# ans = optimr output -- minor modification
# par --    The best set of parameters found.
# value	--   The value of fn corresponding to par.
# counts	-- A two-element integer vector giving the number of calls to fn and gr 
#            respectively. This excludes those calls needed to compute the Hessian
#             even though the opm() result will have these counts
# convergence	-- An integer code. 0 indicates successful completion
# message	-- A character string which for optim() or optimr() may give additional 
#            information returned by the optimizer, or NULL.
#            Here will be "Result of conversion from opm() result"
# 
## want opmmat <- matrix(NA, nrow=nmeth, ncol=npar+8)

optimr2opm <- function(ans, opmmat, xtime=NA){
# ans is the optimr structure solution
# opmmat is a matrix form of the opm() output object (NOT the summary() result)
   npar<-length(ans$par)
   pstring<-NULL
   if (is.null(pstring)) {
      for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
   } 
   cnames <- c(pstring, "value", "fevals", "gevals", "hevals", "convergence", "kkt1", "kkt2", "xtime")
   print(cnames)
   kkt1<-NA
   kkt2<-NA # could add these later
   fevals<-optsp$kfn
   gevals<-optsp$kgr
   hevals<-optsp$khe # NOTE: hope these have been updated
#   cat("fevals = ",fevals," counts:"); print(ans$counts)
   # Could check
#   if ( fevals != ans$counts[[1]] || gevals != ans$counts[[2]]) stop("optimr2opm: Counts mismatch")
   addvec <- c(ans$par, ans$value, fevals, gevals, hevals, ans$convergence, kkt1, kkt2, xtime)
   names(addvec)<-cnames
   statusvec <- attr(ans$par, "status")
   if (!exists("opmmat")) {
      npopm<-attr(opmmat, "npar")
      msg<-paste("optimr2opm: parameter vector length missmatch: optimr->",npar," opm->",npopm,sep='')
      if (npar != npopm) stop(msg)
      opmmat <- matrix(addvec, ncol = length(addvec))
#      colnames(opmmat)<-cnames
      statusmat <- matrix(statusvec, ncol=npar)
   } else
   {
      statusmat <- attr(opmmat, "statusmat")
      opmmat <- rbind(opmmat, addvec)
      statusmat <- rbind(statusmat, statusvec)
   }
   # print(opmmat)
   attr(opmmat, "statusmat")<-statusmat
   opmmat
}  