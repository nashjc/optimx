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
  addvec <- c(ans$par, ans$value, fevals, gevals, hevals, 
              ans$convergence, kkt1, kkt2, xtime)
  names(addvec)<-cnames
  statusvec <- attr(ans$par, "status")
  havemat<-exists("opmmat")
  methn<-attr(ans$value, "method")
  if (! havemat) {
    opmmat <- matrix(addvec, ncol = length(addvec))
    colnames(opmmat)<-cnames
    statusmat <- matrix(statusvec, ncol=npar)
    row.names(opmmat)[1]<-methn
    nrow<-1
  } else
  { 
    npopm<-attr(opmmat, "npar")
    msg<-paste("optimr2opm: parameter vector length missmatch: optimr->",npar," opm->",npopm,sep='')
    if (npar != npopm) stop(msg)
    statusmat <- attr(opmmat, "statusmat")
    opmmat <- rbind(opmmat, addvec)
    nrow <- dim(opmmat)[1] # use new dimentions
    row.names(opmmat)[ nrow ] <- methn
    statusmat <- rbind(statusmat, statusvec)
  }
  kktres <- list(gmax=NA, evratio = NA, kkt1=NA, kkt2=NA, 
                 hev=rep(NA,npar), ngatend=rep(NA, npar), 
                 nhatend=rep(NA, npar*npar))
  # put together results
  ans$xtimes <- xtime # just in case
  amax<-attr(ans,"maximize")
  if (is.null(ans$message)) {amsg<-"none"} else {amsg<-ans$message}
  if (havemat) {
    odetails<-attr(opmmat, "details")
    omax<-attr(opmmat,"maximize")
    odetails<-rbind(odetails, list(method=methn, ngatend=as.numeric(kktres$ngatend), 
                       nhatend=as.numeric(kktres$nhatend), hev=as.numeric(kktres$hev), 
                       message=amsg))
  } else
  { # no matrix yet  
     odetails <- list(method=methn, ngatend=kktres$ngatend, 
                    nhatend=kktres$nhatend, hev=kktres$hev, message=amsg)
     cat("str(odetails):\n")
     print(str(odetails))
  }  
  row.names(odetails)[ nrow ] <- methn
  attr(opmmat, "details") <- odetails
  attr(opmmat, "status")<-statusmat
  answer <- structure(opmmat, details = odetails, maximize = NA,
                      npar = npar, class = c("opm", "data.frame"))
  attr(answer,"fname") <- "optimr2opm"
  attr(answer,"maximize")<-amax
  answer
} 