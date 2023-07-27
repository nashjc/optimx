opm <- function(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=FALSE,
            control=list(),
             ...) {
  npar <- length(par)
  pstring<-names(par) # the names of the parameters
  npar <- length(par)
  ctrl <- ctrldefault(npar)
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
  }
  control <- ctrl # note the copy back! control now has a FULL set of values
  ## 180706: Should we try to streamline?
  if(control$trace > 0) cat("opm: wrapper to call optimr to run multiple optimizers\n")

  fnscale <- 1.0 # default to ensure defined and MINIMIZING
  if (! is.null(control$maximize)){ 
      if ( control$maximize ) {fnscale <- -1.0} 
  }
  else { # control$maximize is NULL, so control$fnscale defines behaviour
      fnscale <- control$fnscale # default is 1.0
      if (fnscale < 0) control$maximize<-TRUE # reset maximize if it was null
      # reset may be needed for kkt check later in opm.
  } # control$maximize has precedence over control$fnscale
  control$fnscale <- fnscale # to ensure set again

  allmeth <- control$allmeth[ - which(control$allmeth %in% control$weakmeth) ]
  # 160628: uobyqa removed as it fails hobbs from 1,1,1 unscaled

  bdmeth <- control$bdmeth

  maskmeth <- control$maskmeth
  # Masks: As at 2016-6-28 did NOT provide for masks in abandoned package optimr

  bmtst <- bmchk(par, lower=lower, upper=upper)
  # Check for onbound and "nmkb" %in% method list
  if (bmtst$onbound && ("nmkb" %in% method)) { # remove nmkb from method
    method <- method[ - which(method == "nmkb")]
    if (control$trace > 0) cat("Start is on bound. Methods include nmkb.\n")
    warning("Removing nmkb from method list in opm() because start is on bound")
  } 
  if (length(method) < 1) stop("No methods requested for opm()")
  control$have.bounds <- bmtst$bounds # and set a control value
  bdmsk <- bmtst$bdmsk # Only need the masks bit from here on
  # These are set free (1) or set -1 for upper bounds, -3 for lower bounds
  # At this stage should NOT have masks (Or could they be added if upper=lower by bmchk
  control$have.masks <- any(bdmsk == 0)

  # 20230614 -- MOST method addition
  if (control$all.methods) method <- allmeth # set list of methods
  if (length(method) == 1) {
      if (method == "ALL") {
          control$all.methods <- TRUE
          method <- allmeth # and change to a vector
      } else {
        if (method == "MOST"){ # has to be inside or method now length > 1
          method <- control$mostmeth # set to most methods
        }
      }
  }
  method <- unique(method) # in case user has duplicates
  if (length(method) > 1) {
       if (control$have.masks) { method <- intersect(method, maskmeth) } # ?? do we want this?
       else { if (control$have.bounds) { method <- intersect(method,bdmeth) }
              else { method <- intersect(method, allmeth) }
       }
       if ( is.null(hess) ) { # remove snewton and snewtonm when no hessian
          if ( "snewton" %in% method ) {
              method <- method[-which(method == "snewton")]
              warning("'snewton' removed from 'method' -- no hess()")
          }
          if ( "snewtonm" %in% method ) {
              method <- method[-which(method == "snewtonm")]
              warning("'snewtonm' removed from 'method' -- no hess()")
          }
       }
  }
  # 20220221: fixup for methods NOT suitable for bounds
  method <- unlist(method) # ?? needed?
  if (control$have.bounds) method <- method[which(method %in% bdmeth)]
  # end fixup
  nmeth <- length(method)
  if (nmeth < 1) stop("No suitable methods requested for opm()")

  if (is.null(pstring)) {
      for (j in 1:npar) {  pstring[[j]]<- paste("p",j,sep='')}
  } 
  cnames <- c(pstring, "value", "fevals", "gevals", "hevals", "convergence", "kkt1", "kkt2", "xtime")
  ans.ret <- matrix(NA, nrow=nmeth, ncol=npar+8) # add hevals 230619
  ans.ret <- data.frame(ans.ret)
  ans.status <- matrix(" ",nrow=nmeth, ncol=npar)
#  ans.ptype <- rep(" ",npar)
  colnames(ans.ret)<-cnames
  row.names(ans.ret)<-method
  ans.details <- list()
  if (control$trace > 2) {
     cat("width of ans.ret =", npar+8,"\n")
     print(dim(ans.ret))
  }
  for (i in 1:nmeth) {
    meth <- method[[i]] # extract the method name. Note double brackets or get list of length 1.
    if (control$trace > 0) cat("Method: ",meth,"\n")
    # Note: not using try() here
    if (is.character(gr) && (control$trace>0)) 
           cat("Using numerical gradient ",gr," for method ", meth,"\n")
    time <- system.time(ans <- optimr(par, fn, gr, hess=hess, method=meth, lower=lower, upper=upper, 
           hessian=hessian, control=control, ...))[1]
    if (control$trace > 2) print(ans)
    # add to list

## --------------------------------------------
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
      if (control$trace>0) { cat("Post processing for method ",meth,"\n") }
#      cat("opm - post processing: ans$convergence=",ans$convergence,"\n")
      if (exists("ans$message")) {
           amsg<-ans$message
           ans$message <- NULL # Safety. Do we need this?
      } else { amsg <- "none" }
      ngatend <- NA
      nhatend <- NA
      hev <- NA
      ans$gevals <- ans$scounts[2]
      ans$fevals <- ans$scounts[1]
      ans$hevals <- ans$scounts[3]
      ans$kkt1<-NA
      ans$kkt2<-NA
      kktres <- list(gmax=NA, evratio = NA, kkt1=NA, kkt2=NA, 
                     hev=rep(NA,npar), ngatend=NA, nhatend=NA)
      if ( control$save.failures || (ans$convergence < 1) ){
           # Save soln if converged or directed to save
          if ((control$trace > 0) && (ans$convergence==0)) cat("Successful convergence! \n") 
          if (control$trace > 3) cat("ans$convergence=",ans$convergence,"\n") # ??
          if ((control$kkt || hessian) && (ans$convergence < 9000)) { 
             # chg 160917 for no gradient, 220330 for snewton solve failure
             wgr <- gr
             if (is.null(wgr)) wgr <- control$defgrapprox
             kktres <- kktchk(ans$par, fn, wgr, hess=NULL, upper=NULL, lower=NULL, 
                    maximize=control$maximize, control=control, ...) 
             ans$kkt1<-as.logical(kktres$kkt1)
             ans$kkt2<-as.logical(kktres$kkt2)
          }
# put together results
          ans$xtimes <- time
          # Do we want more information saved?
          if (control$trace > 1) { 
		cat("Save results from method ",meth,"\n") 
	  	print(ans)
	  }
	  if (control$trace > 2) { 
             cat("Assemble the answers\n") 
             cat("ans.ret now\n")
             print(ans.ret)
          }
          addvec <- c(ans$par, ans$value, ans$fevals, ans$gevals, ans$hevals,
                              ans$convergence, ans$kkt1, ans$kkt2, ans$xtimes)
	  if (control$trace > 2) { 
             cat("length addvec = ", length(addvec),"\n")
             print(addvec)
          }
          ans.ret[i, ] <- addvec
#          ans.ptype[i]  <- attr(ans$value,"ptype")
#          cat("ans.ptype:\n"); print(ans.ptype)
          statusvec <- attr(ans$par, "status")
#          cat("statusvec:");print(statusvec)
          ans.status[i, ] <- statusvec
      }  ## end post-processing of successful solution
      ans.details<-rbind(ans.details, list(method=meth, ngatend=kktres$ngatend, 
             nhatend=kktres$nhatend, hev=kktres$hev, message=amsg))
      # 1303234 try list() not c()
      row.names(ans.details)[[i]]<-meth
    } # End loop  ## end loop over method (index i)
    ansout <- NULL # default if no answers
    if (length(ans$par) > 0) { # cannot save if no answers
	ansout <- ans.ret # Don't seem to need drop=FALSE
        attr(ansout, "details")<-ans.details
        attr(ansout, "status")<-ans.status
#        attr(ansout, "ptype")<-ans.ptype
        ansout[, "kkt1"] <- as.logical(ansout[, "kkt1"])
        ansout[, "kkt2"] <- as.logical(ansout[, "kkt2"])
    }
    ansout # return(ansout)
    answer <- structure(ansout, details = ans.details, maximize = control$maximize,
            npar = npar, class = c("opm", "data.frame"))

} ## end of opm

