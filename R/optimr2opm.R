optimr2opm <- function(ans, opmmat){
# ans is an optimr structure solution
# opmmat is a matrix form of the opm() output object (NOT the summary() result)
# This will be created if it doesn't exist.
   npar<-length(ans$par)
   pstring<-NULL
   if (is.null(pstring)) {
      for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
   } 
   cnames <- c(pstring, "value", "fevals", "gevals", "hevals", "convergence", "kkt1", "kkt2", "xtime")
   kkt1<-NA
   kkt2<-NA # could add these later
   fevals<-optsp$kfn
   gevals<-optsp$kgr
   hevals<-optsp$khe # NOTE: hope these have been updated
   if (is.null(ans$xtime)) {xtime <- NA} else { xtime <- ans$xtime }
   addvec <- c(ans$par, ans$value, fevals, gevals, hevals, ans$convergence, kkt1, kkt2, xtime)
   names(addvec)<-cnames
   statusvec <- attr(ans$par, "status")
   if (!exists("opmmat")) {
      opmmat <- matrix(addvec, ncol = length(addvec))
      colnames(opmmat)<-cnames
      statusmat <- matrix(statusvec, ncol=npar)
      row.names(opmmat)[1]<-attr(ans$value, "method")
   } else
   { npopm<-attr(opmmat, "npar")
     msg<-paste("optimr2opm: parameter vector length missmatch: optimr->",npar," opm->",npopm,sep='')
     if (npar != npopm) stop(msg)
     statusmat <- attr(opmmat, "statusmat")
      opmmat <- rbind(opmmat, addvec)
      row.names(opmmat)[dim(opmmat)[1] ] <- attr(ans$value, "method")
      statusmat <- rbind(statusmat, statusvec)
   }
   kktres <- list(gmax=NA, evratio = NA, kkt1=NA, kkt2=NA, 
                  hev=rep(NA,npar), ngatend=NA, nhatend=NA)
   # put together results
   ans$xtimes <- xtime # just in case
   ans.details<-rbind(ans.details, list(method=meth, ngatend=kktres$ngatend, 
             nhatend=kktres$nhatend, hev=kktres$hev, message=amsg))
   row.names(ans.details)[[i]]<-meth
   # print(opmmat)
   attr(opmmat, "statusmat")<-statusmat
   opmmat
} 