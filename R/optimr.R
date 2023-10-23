optimr <- function(par, fn, gr=NULL, hess=NULL, method=NULL, lower=-Inf, upper=Inf, 
            hessian=FALSE, control=list(), ...) {

  fname <- as.list(sys.call())$fn # FAILED rlang::as_name(as.list(sys.call())$fn)
  if (length(method) > 1) stop("optimr requires single method")
  if (is.null(method)) method <- control$defmethod # Set a default method
  fn1 <- function(par) fn(par,...) # avoid dotarg clashes
  gr1 <- if (!is.null(gr) && !is.character(gr)) function(par) gr(par,...) # ?? will this fail for quoted names
  he1 <- if (!is.null(hess) && !is.character(hess)) function(par) hess(par,...) 

  npar <- length(par)
  ctrl <- ctrldefault(npar)
  ncontrol <- names(control)
  nctrl <- names(ctrl)
  sctrl <- list() # null list of special controls. Only for single methods!
  for (onename in ncontrol) {
     if (onename %in% nctrl) {
       if (! is.null(control[onename]) || ! is.na(control[onename]) )
       ctrl[onename]<-control[onename]
     }
     else { # these are special controls
       sctrl[onename] <- control[onename]
     }
  }
  mcontrol <- list() # define the control list
  nsctrl<-length(sctrl)
  if (nsctrl > 0) {
     warning("Special controls present for optimr with method ",method)
     if (! is.null(control$trace) && (control$trace > 0)) { 
        cat("Special controls present for optimr with method ",method,":")
        for (ic in 1:nsctrl){
            namt<-names(sctrl)[ic]; vt <- unlist(sctrl[ic]); cat(namt," = ",vt,"\n")
        }
     }
     mcontrol <- sctrl
  }
  if (is.null(fname)) fname <- "(no_name)"
  control <- ctrl # note the copy back! control now has a FULL set of values
  # select numerical approx for jacobian and hessian
  jacf <- getFromNamespace("jacobian", control$jacpkg)
  hesf <- getFromNamespace("hessian", control$hesspkg)
  if (control$hesspkg != "numDeriv") warning("Using Hessian approximation from ",control$hesspkg)
  if (control$jacpkg != "numDeriv") warning("Using Jacobian approximation from ",control$hesspkg)
  
  outmethod <- checksolver(method, control$allmeth, control$allpkg) # there will only be one! 
  if (is.null(outmethod)) {
		if (control$trace > 0) cat("Solver ",method," missing\n")
		ans<-list() # ans not yet defined, so set as list
    ans$convergence <- 8888 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
    attr(ans$par, "status")<-rep("?",npar)
	  ans$counts[1] <- NA # save function and gradient count information
	  ans$counts[2] <- NA # save function and gradient count information
	  ans$message <- paste("Missing method ",method)
    ans$hessian <- NULL
    return(ans) # can't proceed without solver
  }

# Check if bounded
  if (is.null(control$trace)) control$trace <- 0 # Fix 220221
# ?? check if method is a bounded method?

  bdmsk <- bmchk(par, lower=lower, upper=upper, shift2bound=TRUE, trace=control$trace)
  if (! bdmsk$feasible) stop("infeasible starting parameters")
#  if (bdmsk$parchanged) {
#      warning("Parameter(s) changed to nearest bounds")
#      par <- bdmsk$bvec
#  }
  control$have.bounds <- bdmsk$bounds # and set a control value
  if (control$have.bounds && !(method %in% control$bdmeth)) stop("Bounded problem with unsuitable method")

  orig.method <- method
  if (!is.null(gr) && !is.character(gr)) { orig.gr <- gr1 } else {orig.gr <- gr }
  orig.fn <- fn1

# ?? put in checkgrad here somehow??

  if (is.null(hessian) ){
     savehess <- FALSE
  } else { savehess <- hessian } # logical -- whether to save hessian 

  if (is.null(control$parscale)) { 
        pscale <- rep(1,npar)
        if(control$trace > 0) { cat("Unit parameter scaling\n") }
  } else { 
        pscale <- control$parscale 
        if(control$trace > 0) {
          cat("Parameter scaling:")
          print(pscale)
        }
  }
  spar <- par/pscale # scaled parameters
  slower <- -Inf
  supper <- Inf # to ensure defined
  if (control$have.bounds) {
    slower <- lower/pscale
    supper <- upper/pscale
  }
  if (is.null(control$fnscale)) {
     fnscale <- 1.0 # default to ensure defined and MINIMIZING
  }
  else {
     if (control$fnscale != 1.0) warning("optimr: use maximize control rather than fnscale")
     fnscale <- control$fnscale
  }
  if (! is.null(control$maximize)){ 
      if ( control$maximize ) {fnscale <- -1.0} 
  }
  else { # control$maximize is NULL, so control$fnscale defines behaviour
      fnscale <- control$fnscale # default is 1.0
      if (fnscale < 0) control$maximize<-TRUE # reset maximize if it was null and maximizing
  } # control$maximize has precedence over control$fnscale
  control$fnscale <- fnscale # to ensure set again

  efn <- function(spar) { # dotargs removed 230619
      optsp$kfn<-optsp$kfn+1 # counter
      # rely on pscale being defined in this enclosing environment
      par <- spar*pscale
      val <- fn1(par) * fnscale
      if (is.complex(val) || ! is.finite(val) ) val <- control$bigval # ?? maybe sqrt(bigval)?
      val
  }

  if (is.character(gr)) { # approximation to gradient
     if ( ! (gr %in% control$grapprox)) stop(gr," is not a valid gradient approximation code")
     if (control$trace>0) cat("Using numerical approximation '",gr,"' to gradient in optimr()\n")
     egr <- function(spar){ 
       optsp$kgr<-optsp$kgr+1 # counter
       if (control$trace > 2) {
         cat("par:"); print(par)
         cat("fnscale =",fnscale,"  pscale="); print(pscale)
         cat("gr:"); print(gr)
       }
       par <- spar*pscale
       result <- do.call(gr, list(par, userfn=fn1)) * pscale * fnscale # dotargs removed 230619
       ## need fnscale.
    }
  } else { 
    if (is.null(gr)) {egr <- NULL}
    else {
       egr <- function(spar) {
         optsp$kgr<-optsp$kgr+1 # counter
         par <- spar*pscale
         result <- gr1(par) * pscale * fnscale
       }
    }
  } # end egr definition

  orig.hess<-NULL # to ensure defined
  if (is.null(hess)) { 
     ehess <- NULL 
  }
  else { # check for approx 
    if (is.character(hess)){
#      cat("hess is character ==",hess,"\n")
      if (hess != "approx") {
          stop("Undefined character hessian -- hess =",hess)
      } else {
        if (is.null(gr) || is.character(gr) ) { # use hesf
          if (control$trace>0) cat("Using numerical approx. ",control$hesspkg,"::hessian in optimr()\n")
          if (control$trace > 2) {
            cat("par:"); print(par)
            cat("fnscale =",fnscale,"  pscale="); print(pscale)
          }
          ehess <- function(spar){
            optsp$khe<-optsp$khe+1 # counter
            par <- spar*pscale
            result <- hesf(fn1, par) * pscale * pscale * fnscale  ## need fnscale.
          } # end ehess for fn but no gr and hess == "approx"
          tgr <- "grcentral" ## Need to add egr, use "grcentral"
          if (control$trace>0) cat("Approx hessian: replace null gradient with ",tgr," in optimr()\n")
          egr <- function(spar){ 
            optsp$kgr<-optsp$kgr+1 # counter
            par <- spar*pscale
            result <- do.call(tgr, list(par, userfn=fn1)) * pscale * fnscale ## need fnscale.
          }
        }
        else {  # have gr, hess=="approx",  so use jacobian of egr
          if (control$trace>0) cat("Using numerical approx. numDeriv::jacobian of gr in optimr()\n")
          if (control$trace > 2) {
            cat("par:"); print(par)
            cat("fnscale =",fnscale,"  pscale="); print(pscale)
          }
          ehess <- function(spar){
            optsp$khe<-optsp$khe+1 # counter
            par <- spar*pscale
            result <- jacf(gr1, par) * pscale * pscale * fnscale  ## need fnscale.
          } # end ehess for fn and gr but no hess
        } # end gr and fn no hess
      }
    } # end character hess
    else { # Not character hess, and not null
       if (is.null(egr)) stop("Null egr but non-null hess -- optimr failure")
       orig.hess<-hess
       ehess <- function(spar) {
         optsp$khe<-optsp$khe+1 # counter
         par <- spar*pscale
         result <- he1(par) * pscale * pscale * fnscale # dotargs removed 230619
         result
       }
    } # end not null hess 
  } # end check for null hess
  
  if (control$trace > 2) {
    cat("spar = "); print(spar)
    cat("efn = ",efn(spar),"\n"); # dotargs removed 230619
    cat("egr:"); if (is.null(egr)) print("is NULL") else print(egr(spar)) # dotargs removed 230619
    cat("ehess:"); if (is.null(ehess)) print("is NULL") else print(ehess(spar)) # dotargs removed 230619 
  }
  
  # Set counters for fn, gr, hess. NOTE: do this AFTER display, or get +1 counts
  optsp$kfn <- 0
  optsp$kgr <- 0
  optsp$khe <- 0

  fghfn <- NULL # ensure fn for nlm is undefined unless we need it
  fgfn <- NULL # ensure fn for nlm is undefined unless we need it
  if ("nlm" == method) {
    fghfn <- function(spar){ # dotargs removed 230619
      f <- efn(spar) # dotargs removed 230619
      if (is.null(egr)) {g <- NULL} else {g <- egr(spar)} # dotargs removed 230619 
      attr(f,"gradient") <- g
      if (is.null(ehess)) { h <- NULL } else {h <- ehess(spar)} # dotargs removed 230619
      attr(f,"hessian") <- h
      f
    }
  } # end definition of fghfn
  if ("Rtnmin" == method) {
    fgfn <- function(spar){ # dotargs removed 230619
      f <- efn(spar) # dotargs removed 230619
      if (is.null(egr)) {g <- NULL} else {g <- egr(spar)} # dotargs removed 230619 
      attr(f,"gradient") <- g
      f
    }
  } # end definition of fgfn

## Masks 
   maskmeth <- control$maskmeth
   msk <- bdmsk$bdmsk # Only need the masks bit from here on
   if (any(msk == 0) ) {
      if ( !(method %in% maskmeth) ) {
         stopmsg <- paste("Method ",method," cannot handle masked (fixed) parameters")
         stop(stopmsg)
      }
      if (control$trace > 0) cat("Masks present\n")
   }

# replacement for optim to minimize using a single method

# time is in opm(), but NOT here
# The structure has   par, value, counts, convergence, message, hessian

# expand bounds
  if (length(lower) == 1 && is.finite(lower) ) lower<-rep(lower,npar)
  if (length(upper) == 1 && is.finite(upper) ) upper<-rep(upper,npar)
# set solver controls
  mcontrol$maximize <- NULL # efn, egr have the control
# Methods from optim()
  if (method == "SANN") warning("optim::SANN is inappropriate for use with optimx package!")
  if (method == "Nelder-Mead" || 
      method == "BFGS" || 
      method == "L-BFGS-B" || 
      method == "CG" || 
      method == "SANN") {
      # Take care of methods from optim(): Nelder-Mead, BFGS, L-BFGS-B, CG
      mcontrol$trace <- control$trace 
      mcontrol$parscale <- NULL # using efn/egr. DO NOT USE fnscale, parscale for optim methods
      mcontrol$fnscale <- NULL
      mcontrol$maxit <- control$maxit 
      if (! is.null(control$maxit)) {mcontrol$maxit <- control$maxit}
# Note: hessian always FALSE in these calls. But savehess may recover it.
#        cat("Before optim() call - control$have.bounds =",control$have.bounds,"\n")
      if(nsctrl > 0) { 
          mcontrol<-c(mcontrol,sctrl) 
          cat("controls added:"); print(sctrl)
          print(mcontrol)
      }
      if (control$have.bounds) {
        if (method != "L-BFGS-B") {
            errmsg <- "optim() can only handle bounds with L-BFGS-B\n"
            if (control$trace > 0) cat(errmsg,"\n")
            ans <- list()
            class(ans)[1] <- "try-error"
            warning("optimr: optim() with bounds ONLY uses L-BFGS-B")
        } else {
            if(control$trace > 3) cat("L-BFGS-B start: ")
            ans <- try(optim(par=par, fn=efn, gr=egr, 
                      lower=lower, upper=upper, method="L-BFGS-B", hessian=FALSE, 
                       control=mcontrol)) # remove , ... 231016
            if(control$trace > 3) cat(" ans$value=",ans$value,"\n")
          }
        } else {
          ans <- try(optim(par=par, fn=efn, gr=egr, 
                method=method, hessian=FALSE, control=mcontrol)) # remove , ... 231016
        }
        if (inherits(ans,"try-error")) { # bad result -- What to do?
		  ans<-list() # ans not yet defined, so set as list
          ans$convergence <- 9999 # failed in run
          errmsg <- "optim method failure\n"
          if (control$stoponerror) stop(errmsg)
          if (method != "L-BFGS-B") 
              errmsg <- paste("optim() with bounds ONLY uses L-BFGS-B: ", errmsg)
		  if (control$trace>0) cat(errmsg)
	      ans$value <- control$badval
		  ans$par<-rep(NA,npar)
	      ans$counts[1] <- NA # save function and gradient count information
	      ans$counts[2] <- NA # save function and gradient count information
	      ans$message <- errmsg
        } # otherwise ans is OK and we return it
        ## return(ans) # to ensure we return
      }   # end if using optim() methods
## --------------------------------------------
      else if (method == "nlminb") {
        # Here we use portLib routine nlminb rather than optim as our minimizer
        mcontrol$iter.max<-mcontrol$maxit # different name for iteration limit in this routine
        mcontrol$maxit<-NULL # and we null it out
        mcontrol$abs.tol <- 0 # To fix issues when minimum is less than 0. 20100711
        mcontrol$eval.max <- control$maxfeval
	if ( is.null(control$trace) || is.na(control$trace) || control$trace == 0) { 
		mcontrol$trace = 0
	} else { 
		mcontrol$trace = 1 # this is EVERY iteration. nlminb trace is freq of reporting.
	}
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        ans <- try(nlminb(start=spar, objective=efn, gradient=egr, hessian=ehess, lower=slower, 
		upper=supper, control=mcontrol))
        if (! inherits(ans, "try-error")) {
		# Translate output to common format and names
        	ans$value<-ans$objective
                ans$par <- ans$par*pscale
	        ans$objective<-NULL
	        ans$counts[1] <- ans$evaluations[1]
        	ans$counts[2] <- ans$evaluations[2]
		ans$evaluations<-NULL # cleanup
	        ans$iterations<-NULL
                ans$hessian <- NULL
	} else { # bad result -- What to do?
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		if (control$trace>0) cat("nlminb failure\n")
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- "nlminb failed" # 180318 change from NULL
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using nlminb
## --------------------------------------------
      else if (method == "nlm") { # Use stats package nlm routine
#        if (is.null(gr)) { stop("optimr -- nlm -- we do not allow gr = NULL") }
	if (! is.null(control$maxit) ) {iterlim <- control$maxit }
        else { iterlim <- 100 }
	print.level <- 0 # ??trace??
        errmsg <- NULL
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (control$have.bounds) {
              if(control$trace > 0) cat("nlm cannot handle bounds\n")
              errmsg <- "nlm cannot handle bounds\n"
            ##  stop("nlm tried with bounds")
            ans <- list()
            class(ans)[1] <- "try-error"
        } else {
          if (! is.null(control$trace) && (control$trace > 0) ) {print.level <- 2 } 
#          ans <- try(nlm(f=fghfn, p=spar, iterlim=iterlim, print.level=print.level, ...))
           ans <- try(nlm(f=fghfn, p=spar, iterlim=iterlim, print.level=print.level))
        }
        if (! inherits(ans, "try-error")) {
		if (ans$code == 1 || ans$code == 2 || ans$code == 3) ans$convergence <- 0
		if (ans$code == 4) ans$convergence <- 1
                if (ans$code == 5) ans$convergence <- 5
        	# Translate output to common format
		ans$value <- ans$minimum
		ans$minimum <- NULL
                ans$par <- ans$estimate*pscale
		ans$estimate <- NULL
                ans$counts[2] <- ans$iterations
                ans$counts[1] <- NA
        	ans$iterations <- NULL
                ans$hessian <- NULL
                ans$gradient <- NULL # We lose information here
                ans$message <- paste("nlm: Convergence indicator (code) = ",ans$code)
                ans$code <- NULL
	} else {
		if (control$trace > 0) cat("nlm failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- "nlm failed" # 180318 change from NULL
                ans$hessian <- NULL
        }
        print.level <- NULL # clean up
        ## return(ans)
      } # end if using nlm
## --------------------------------------------
      else if (method == "ncg") { # Use ncg routine 
        mcontrol$trace <- control$trace
        mcontrol$maxit <- control$maxit # 151217 JN
        mcontrol$stepredn <- control$stepredn # 220217
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (! is.null(egr)) { #?? why not slower, supper?
   	     ans <- try(ncg(par=spar, fn=efn, gr=egr, bds=bdmsk, control=mcontrol))
        }
        if (!is.null(egr) && !inherits(ans, "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (control$trace>0) {
                    cat("ncg failed for current problem \n")
                    if(is.null(egr)) cat("Note: ncg needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                if(is.null(egr)) {
                   ans$message <- "Must specify gradient function for ncg"       
                   ans$convergence <- 9998 # for no gradient where needed
                }
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using ncg
## --------------------------------------------
      else if (method == "Rcgmin") { # Use Rcgmin routine (ignoring masks)
        mcontrol$trace <- control$trace
        mcontrol$maxit <- control$maxit # 151217 JN
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (! is.null(egr)) {
  	  if (control$have.bounds) { # 151220 -- this was not defined
            # 170919 -- explicit reference to package
   	    ans <- try(Rcgminb(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol, ...))
	  } else {
#   	     ans <- try(Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol, ...))
   	     ans <- try(Rcgminu(par=spar, fn=efn, gr=egr, control=mcontrol))
	  }
        }
        if (!is.null(egr) && !inherits(ans, "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (control$trace>0) {
                    cat("Rcgmin failed for current problem \n")
                    if(is.null(egr)) cat("Note: Rcgmin needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                if(is.null(egr)) {
                   ans$message <- "Must specify gradient function for Rcgmin"       
                   ans$convergence <- 9998 # for no gradient where needed
                }
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rcgmin
## --------------------------------------------
      else if (method == "nvm") { # Use nvm routine 
        mcontrol$trace <- control$trace
        mcontrol$maxit <- control$maxit # 151217 JN
        mcontrol$stepredn <- control$stepredn # 220217
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (! is.null(egr)) { #?? why not slower, supper?
# 230622 update to simplify call
   	     ans <- try(nvm(par=spar, fn=efn, gr=egr, bds=bdmsk, control=mcontrol))
#   	     ans <- try(nvm(par=spar, fn=efn, gr=egr, lower=slower, upper=supper,
#                         bdmsk=bdmsk$bdmsk, control=mcontrol))
#                        bdmsk=bdmsk$bdmsk, control=mcontrol, ...))
        }
        if (!is.null(egr) && !inherits(ans, "try-error")) {
                ans$par <- ans$par*pscale
	        ans$message <- NA        
                ans$hessian <- NULL
                ans$bdmsk <- NULL # clear this
        } else {
		if (control$trace>0) {
                    cat("nvm failed for current problem \n")
                    if(is.null(egr)) cat("Note: nvm needs gradient function specified\n")
                }
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                if(is.null(egr)) {
                   ans$message <- "Must specify gradient function for nvm"       
                   ans$convergence <- 9998 # for no gradient where needed
                }
                ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using nvm
## --------------------------------------------
      else if (method == "Rvmmin") { # Use Rvmmin routine
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
        mcontrol$parscale <- NULL # using efn/egr. DO NOT USE fnscale, parscale
        mcontrol$fnscale <- NULL
	mcontrol$trace <- control$trace # 140902 Note no check on validity of values
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
	if (! is.null(egr)) {
          ans <- try(Rvmmin(par=spar, fn=efn, gr=egr, lower=slower,
                upper=supper, bdmsk=msk, control=mcontrol))
#                upper=supper, bdmsk=msk, control=mcontrol, ...))
        }
        if (control$trace > 2) {
            cat("Rvmmin ans:")
            print(ans)
        }
        if (! is.null(egr) && !inherits(ans, "try-error")) {
            ans$par <- ans$par*pscale
            ans$bdmsk <- NULL
        } else {
            if (control$trace>0) {
                cat("Rvmmin failed for current problem \n")
                if(is.null(egr)) cat("Note: Rvmmin needs gradient function specified\n")
            }
	    ans<-list() # ans not yet defined, so set as list

	    ans$value <- control$badval
	    ans$par<-rep(NA,npar)
	    ans$counts[1] <- NA # save function and gradient count information
	    ans$counts[2] <- NA # save function and gradient count information
	    ans$message <- NULL        
            if(is.null(egr)) {
               ans$message <- "Must specify gradient function for Rvmmin"
               ans$convergence <- 9998 # for no gradient where needed
            }
            ans$hessian <- NULL
        }
        ## return(ans)
      }  ## end if using Rvmmin
## --------------------------------------------
      else if (method == "snewton") { # Use snewton routine
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval # changed from maxfevals 180321
       	mcontrol$trace <- control$trace # 140902 Note no check on validity of values
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
       	ans<-list(par=NA, value=NA, counts=NA, convergence=-1,
                  message=NA, hessian=NULL) # ans not yet defined, so set as list
       	if ( (! is.null(egr)) && (! is.null(ehess)) ) {
          if (control$have.bounds) { # 170919 make package explicit
             stop("snewton does not handle bounds") 
          } 
#       	  tans <- try( snewton(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol,...))
        	  tans <- try( snewton(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol))
          if (control$trace>1) {
              cat("snewton returns tans:")
              print(tans)
          }
          if  (inherits(tans, "try-error")) { 
             ans$message <- "snewton failed"
             ans$convergence <- 9999
             if (control$trace>0) {
               cat(ans$message,"\n")
             }
          }
       	} else {
       	   if(is.null(egr)) {
       	     ans$message <- "Must specify gradient function for snewton"
       	     ans$convergence <- 9998 # for no gradient where needed
       	     warning("Note: snewton needs gradient function specified")
       	   }
       	   if(is.null(ehess)) {
       	     ans$message <- "Must specify Hessian function (hess) for snewton"
       	     ans$convergence <- 9997 # for no Hessian where needed
       	     warning("Note: snewton needs Hessian function (hess) specified")
       	   }
       	}
       	if (ans$convergence > 9996){
       	       ans$value <- control$badval
               ans$par<-rep(NA,npar)
               ans$counts[1] <- NA # save function and gradient count information
               ans$counts[2] <- NA # save function and gradient count information
               # Note: in optim() no provision for hessian count
               ans$hessian <- NULL
               if (control$trace>1) {
                  cat("snewton falure ans:")
                  print(ans)
               }
        } else { # have an answer
              ans$par <- tans$par*pscale
              ans$value <- tans$value
              attr(ans, "gradient") <- tans$grad
              if(hessian) ans$hessian <- tans$Hess
              ans$counts[1] <- tans$counts$nfn
              ans$counts[2] <-  tans$counts$ngr
              ans$message <- tans$message 
              ans$convergence <- tans$convcode
              tans <- NULL # probably unnecessary, but for safety
              if (control$trace>1) {
                   cat("rejigged ans:")
                   print(ans)
              }
             } # end have answer
         ## return(ans)
      }  ## end if using snewton
  ## --------------------------------------------
  else if ((method == "snewtm") || (method=="snewtonm")) { # Use snewtm/snewtonm routine
    mcontrol$maxit <- control$maxit
    mcontrol$maxfeval <- control$maxfeval # changed from maxfevals 180321
    mcontrol$trace <- control$trace # 140902 Note no check on validity of values
    ans<-list(par=NA, value=NA, counts=NA, convergence=-1,
                  message=NA, hessian=NULL) # ans not yet defined, so set as list
    ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
    if ( (! is.null(egr)) && (! is.null(ehess)) ) {
      if (control$have.bounds) { # 170919 make package explicit
        tans <- try( snewtm(par=spar, fn=efn, gr=egr, hess=ehess, bds=bdmsk, control=mcontrol))
#        tans <- try( snewtonm(par=spar, fn=efn, gr=egr, hess=ehess, lower=slower,
#                upper=supper, control=mcontrol))
#                upper=supper, control=mcontrol, ...))
        # added 20220210. masks NOT checked yet??
      }
#      else tans <- try( snewtonm(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol,...))
#       else tans <- try( snewtonm(par=spar, fn=efn, gr=egr, hess=ehess, control=mcontrol))
        else tans <- try( snewtm(par=spar, fn=efn, gr=egr, hess=ehess, bds=bdmsk, control=mcontrol))
      if (control$trace>0) {
           cat("snewtonm returns tans:")
           print(tans)
      }
      if  (inherits(tans, "try-error")) { 
        ans$message <- "snewtonm failed"
        ans$convergence <- 9999
        if (control$trace>0) {
          cat(ans$message,"\n")
        }
      }
    } else {
      if(is.null(egr)) {
        ans$message <- "Must specify gradient function for snewtonm"
        ans$convergence <- 9998 # for no gradient where needed
        warning("Note: snewtonm needs gradient function specified")
      }
      if(is.null(ehess)) {
        ans$message <- "Must specify Hessian function (hess) for snewtonm"
        ans$convergence <- 9997 # for no Hessian where needed
        warning("Note: snewtonm needs Hessian function specified")
      }
    } # end of fails
    if (ans$convergence > 9996){ # Bad solution
       ans$value <- control$badval
       ans$par<-rep(NA,npar)
       ans$counts[1] <- NA # save function and gradient count information
       ans$counts[2] <- NA # save function and gradient count information
       # Note: in optim() no provision for hessian count
       ans$hessian <- NULL
       if (control$trace>1) {
          cat("snewtonm falure ans:")
          print(ans)
       }
    } else { # have an answer
#      cat("copy answer with tans$par:\n") 
#      print(tans$par)              
      ans$par <- tans$par*pscale
      ans$value <- tans$value
      attr(ans, "gradient") <- tans$grad
      if(hessian) ans$hessian <- tans$Hess
      ans$counts[1] <- tans$counts$nfn
      ans$counts[2] <-  tans$counts$ngr
      ans$message <- tans$message 
      ans$convergence <- tans$convcode
      tans <- NULL # probably unnecessary, but for safety
      if (control$trace>1) {
         cat("rejigged ans:")
         print(ans)
      }
    } # end have answer
    ## return(ans)
  }  ## end if using snewtm
  ## --------------------------------------------
      else if (method == "hjn") {# Use JN Hooke and Jeeves
        if (control$trace > 1) { 
           cat("hjn:control$have.bounds =",control$have.bounds,"\n")
           cat("optimr - hjn - msk:"); print(msk)
        }
        mcontrol <- control # copy
        mcontrol$maximize <- NULL # 180327 Cannot maximize within hjn itself.
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        ans <- try(hjn(spar, efn, lower=slower, upper=supper, bdmsk=msk, 
                        control=mcontrol))
#                        control=mcontrol, ...))
        if (! inherits(ans, "try-error")) {
            ans$par <- ans$par*pscale
            ans$message <- "hjn success" 
         } else {
            if (control$trace > 0) cat("hjn failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
            ans$hessian <- NULL
            ans$message <- NA
         }
         ## return(ans)
      }  ## end if using hjn
## --------------------------------------------
      else if (method == "spg") { # Use BB package routine spg as minimizer
        mcontrol$maximize <- NULL # Use external maximization approach
        mcontrol$maxit <- control$maxit
        mcontrol$maxfeval <- control$maxfeval
        if (control$trace > 0) { 
            mcontrol$trace <- TRUE
            if (control$trace > 1) mcontrol$triter <- 1 # default is 10
        } else { mcontrol$trace <- FALSE }
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        ans <- try(BB::spg(par=spar, fn=efn, gr=egr, lower=slower, upper=supper,  
		control=mcontrol))
#		control=mcontrol, ...))
        if (! inherits(ans, "try-error")) {
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval<-NULL # to erase conflicting name
           ans$counts[2] <- ans$iter
           ans$fn.reduction <- NULL # so it does not interfere
           ans$iter<-NULL
           ans$gradient<-NULL # loss of information
        } else { # spg failed
		if (control$trace > 0) cat("spg failed for this problem\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL
                ans$hessian <- NULL
        }
        ## return(ans)
      }  # end if using spg
## --------------------------------------------
      else if (method == "ucminf") {
        ## Use ucminf routine
        if (is.null(control$maxit)) { mcontrol$maxeval <- 500 }  # ensure there is a default value
        else { mcontrol$maxeval <- control$maxit}
        mcontrol$maxit <- NULL # 150427 ensure nulled for ucminf
        errmsg <- NULL
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (control$have.bounds) {
              if (control$trace > 0) cat("ucminf cannot handle bounds\n")
              errmsg <- "ucminf cannot handle bounds\n"
              stop(errmsg)
              ans <- list()
              class(ans)[1] <- "try-error"
          } else {
              uhessian <- 0 # Ensure hessian NOT computed
              ans <- try(ucminf::ucminf(par=spar, fn=efn, gr=egr, 
                   hessian = uhessian,  control=mcontrol))
 #                   hessian = uhessian,  control=mcontrol, ...))
          }
          if (! inherits(ans, "try-error")) {
# From ucminf documentation:  
# convergence = 1 Stopped by small gradient (grtol).
#               2 Stopped by small step (xtol).
#               3 Stopped by function evaluation limit (maxeval).
#               4 Stopped by zero step from line search
#              -2 Computation did not start: length(par) = 0.
#              -4 Computation did not start: stepmax is too small.
#              -5 Computation did not start: grtol or xtol <= 0.
#              -6 Computation did not start: maxeval <= 0.
#              -7 Computation did not start: given Hessian not pos. definite.
#  message: String with reason of termination.
		        if (ans$convergence == 1 
		            || ans$convergence == 2 
		            || ans$convergence == 4) {
         		        ans$convergence <- 0
		        } 
            ans$par <- ans$par*pscale
        	  ans$counts[1] <- ans$info[4]
        	  ans$counts[2] <- ans$info[4] # calls fn and gr together
        	  ans$info <- NULL # to erase conflicting name
        	  ans$nitns <- NULL
            ans$hessian <- NULL
            ans$invhessian.lt <- NULL
		        if (control$trace > 0) cat("ucminf message:",ans$message,"\n")
            } else { # ucminf failed
      		if (control$trace > 0) cat("ucminf failed for this problem\n")
	        ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
	        ans$value <- control$badval
	        ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- errmsg
                ans$hessian <- NULL
            }
            uhessian <- NULL
          ## return(ans)
      }  ## end if using ucminf
## --------------------------------------------
      else if (method == "Rtnmin") { # Use Rtnmin routines 
	if (control$trace>0) {mcontrol$trace <- TRUE } else {mcontrol$trace <- FALSE}
	ans<-list() # ans not yet defined, so set as list
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (is.null(egr)) { ## fixed msg below (referred to lbfgs) 170214
            if (control$trace > 0) cat("Rtnmin MUST have gradient provided\n")
            errmsg <- "Rtnmin MUST have gradient provided"
            class(ans)[1] <- "try-error"            
        } else {
           if (control$have.bounds) {
   	      ans <- try(tnbc(x=spar, fgfun=fgfn, lower=slower,
                   upper=supper, trace=mcontrol$trace))
           } else {
   	      ans <- try(tn(x=spar, fgfun=fgfn, trace=mcontrol$trace))
	   }
        }
        if (inherits(ans,"try-error")) {
           if (control$trace>0) cat("Rtnmin failed for current problem \n")
           ans$convergence <- 9999 # failed in run
	   ans$message <- "Rtnmin failed fo current problem"        
           if (is.null(egr)) {
              ans$convergence <- 9998
              ans$message <- errmsg
              ans$value <- 1234567E20
           } 
	   ans$value <- control$badval
	   ans$par<-rep(NA,npar)
	   ans$counts[1] <- NA # save function and gradient count information
	   ans$counts[2] <- NA 
           ans$hessian <- NULL
        } else {
           ans$par <- ans$xstar*pscale
           ans$xstar <- NULL
           ans$value <- as.numeric(ans$f)
           ans$f <- NULL
           ans$g <- NULL
	   ans$convergence <- ans$ierror
           ans$ierror <- NULL
	   ans$counts[1] <- ans$nfngr
	   ans$counts[2] <- ans$nfngr
           ans$nfngr <- NULL
           ans$hessian <- NULL
	   ans$message <- NA
        }
        ## return(ans)
      }  ## end if using Rtnmin
## --------------------------------------------
      else if (method == "bobyqa") {# Use bobyqa routine from minqa package
        # cat("Entering bobyqa, mcontrol="); print(mcontrol)
       	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        mcontrol$trace <- NULL
        if (is.null(mcontrol$rhobeg)) { # minqa bobyqa defaults 
             mcontrol$rhoend <- NULL # to ensure defaults
        } else {
          if (mcontrol$rhobeg < 0) {  # minqa bobyqa defaults 
             mcontrol$rhobeg <- NULL; mcontrol$rhoend <- NULL # force defaults
          } else {
            if (mcontrol$rhobeg == 0) { # optrimr defaults 
              myrhobeg <- min(supper - slower)/3 # JN 160107 (3), 160125 (5)
              if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
              mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
              mcontrol$rhoend <- 1e-6 * myrhobeg # use 1e-6 scaling as per manual for bobyqa
            }
          } # end else mcontrol$rhobeg >= 0
          # otherwise using avaluable controls
        }
        if (control$trace > 0) {cat("mcontrol for bobyqa:\n"); print(str(mcontrol))}
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        ans <- try(minqa::bobyqa(par=spar, fn=efn, lower=slower,
                upper=supper, control=mcontrol))
        if (! inherits(ans, "try-error")) {
		ans$convergence <- 0
#                if (ans$feval > mcontrol$maxfun) {
#			ans$convergence <- 1 # too many evaluations
#                }
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
                ans$feval <- NULL
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using bobyqa
## --------------------------------------------
      else if (method == "uobyqa") {# Use uobyqa routine from minqa package
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (control$have.bounds) {
            warning("Cannot use uobyqa with bounds")
		if (control$trace > 0) cat("Cannot use uobyqa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
                ## return(ans)
        }
        else {
#          ans <- try(minqa::uobyqa(par=spar, fn=efn, control=mcontrol,...))
           ans <- try(minqa::uobyqa(par=spar, fn=efn, control=mcontrol))
        }
        if (! inherits(ans, "try-error")) {
		ans$convergence <- 0
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
	        ans$counts[2] <- NA
                ans$feval <- NULL
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("uobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using uobyqa
## --------------------------------------------
      else if (method == "newuoa") {# Use newuoa routine from minqa package
        if (control$trace > 1) cat("Trying newuoa\n")
	mcontrol$maxfun <- control$maxfeval
        mcontrol$iprint <- control$trace
        myrhobeg <- min(abs(spar)) # JN 160107 (3), 160125 (5)
        if ((myrhobeg < 1e-8) || ! is.finite(myrhobeg) ) myrhobeg <- 0.5
        mcontrol$rhobeg <- myrhobeg # to avoid 0 when parameters 0
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (control$have.bounds) {
            warning("Cannot use newuoa with bounds")
		if (control$trace > 0) cat("Cannot use newuoa with bounds\n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
                ## return(ans)
        }
#        ans <- try(minqa::newuoa(par=spar, fn=efn, control=mcontrol,...))
        ans <- try(minqa::newuoa(par=spar, fn=efn, control=mcontrol))
        if (! inherits(ans, "try-error")) {
		ans$convergence <- 0
#                if (ans$feval > mcontrol$maxfun) {
#			ans$convergence <- 1 # too many evaluations
#                }
                ans$convergence <- ans$ierr
                ans$ierr <- NULL
                ans$message <- ans$msg
                ans$msg <- NULL
	        ans$counts[1] <- ans$feval
                ans$feval <- NULL
	        ans$counts[2] <- NA
		ans$value<-ans$fval 
                ans$par <- ans$par*pscale
	      	ans$fval <- NULL # not used
                ans$hessian <- NULL
        } else {
		if (control$trace > 0) cat("bobyqa failed for current problem \n")
		ans<-list() # ans not yet defined, so set as list
                ans$convergence <- 9999 # failed in run
		ans$value <- control$badval
		ans$par<-rep(NA,npar)
	        ans$counts[1] <- NA # save function and gradient count information
	        ans$counts[2] <- NA # save function and gradient count information
	        ans$message <- NULL        
                ans$hessian <- NULL
        }
        ans <- unclass(ans) # because minqa does strange things!
        ## return(ans)
      }  ## end if using newuoa
## --------------------------------------------
      else if (method == "nmkb") {# Use nmkb routine from dfoptim package
        if (any(par == lower) || any(par==upper)) {
           if (control$trace>0) cat("nmkb cannot start if on any bound \n")
           warning("nmkb() cannot be started if any parameter on a bound")
           ans <- list() # ans not yet defined, so set as list
           ans$value <- control$badval
           ans$par <- rep(NA,npar)
           ans$convergence <- 9999 # failed in run - !! consider special code for nmkb on bounds
           ans$fevals <- NA 
           ans$gevals <- NA 
           ans$nitns <- NA
           ans$hessian <- NULL
           ans$message <- "nmkb cannot start on bound"
        } else { # ok to proceed with nmkb()
#        if (is.null(control$maxit)) { 
#	   mcontrol$maxfeval <- 5000*round(sqrt(npar+1)) # default at 100215, Change?
#           cat("nmkb maxfeval =",mcontrol$maxfeval
#	}
#        else { 
         mcontrol$maxfeval <- control$maxfeval 
         if (control$trace > 0) { mcontrol$trace <- TRUE } # logical needed, not integer         
         else { mcontrol$trace<-FALSE }
         ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
         if (control$have.bounds) {
            ans <- try(dfoptim::nmkb(par=spar, fn=efn, lower = slower, 
              upper = supper, control=mcontrol))
#              upper = supper, control=mcontrol, ...))
         } else {# 170919 explicit package in call
#            ans <- try(dfoptim::nmk(par=spar, fn=efn, control=mcontrol, ...))
             ans <- try(dfoptim::nmk(par=spar, fn=efn, control=mcontrol))
         }
         if (control$trace > 1) {
            cat("Outputting ans for nmkb:\n")
            print(ans)
         }

        if (! inherits(ans, "try-error")) {
           ans$value <- as.numeric(ans$value)
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval <- NULL
           ans$counts[2] <- NA
      	   ans$nitns <- NA # not used
           # What about 'restarts' and 'message'?!!
           ans$message <- paste(ans$msg," Restarts for stagnation =",ans$restarts)
           ans$restarts <- NULL
           ans$hessian <- NULL
         } else {
           if (control$trace>0) cat("nmkb failed for current problem \n")
           ans <- list(fevals=NA) # ans not yet defined, so set as list
           ans$value <- control$badval
           ans$par <- rep(NA,npar)
           ans$counts[1] <- NA
           ans$counts[2] <- NA
           ans$convergence <- 9999 # failed in run
           ans$message<-"Failed"
           ans$hessian <- NULL
         }
       } # end of check for parameter on bound
       ## return(ans)
     }  ## end if using nmkb
## --------------------------------------------
      else if (method == "hjkb") {# Use hjkb routine from dfoptim package
         if (control$trace > 0) {
            mcontrol$info <- TRUE # logical needed, not integer         
         } else { mcontrol$info <- FALSE }
         mcontrol$maxfeval <- control$maxfeval
         ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
         if (control$have.bounds) {
            ans <- try(dfoptim::hjkb(par=spar, fn=efn, lower = slower, 
                upper = supper, control=mcontrol))
#                upper = supper, control=mcontrol, ...))
         } else {
            ans <- try(dfoptim::hjk(par=spar, fn=efn, control=mcontrol))
         }
         if (! inherits(ans, "try-error")) {
           ans$value <- as.numeric(ans$value)
           ans$par <- ans$par*pscale
           ans$counts[1] <- ans$feval
           ans$feval <- NULL
           ans$counts[2] <- NA
      	   ans$nitns <- NULL # not used
           ans$restarts <- NULL
           ans$hessian <- NULL
           ans$nitns <- NULL # loss of information
         } else {
            if (control$trace>0) cat("hjkb failed for current problem \n")
            ans <- list(value=control$badval, par=rep(NA,npar), message="Failed",
                convergence=9999)
            ans$counts[1]<- NA
            ans$counts[2]<- NA 
            ans$hessian <- NULL
         }
         ## return(ans)
      }  ## end if using hjkb
## --------------------------------------------
      else if (method == "lbfgsb3c") {# Use 2011 L-BFGS-B wrapper
        if (control$trace > 1) cat("lbfgsb3c\n")
        mcontrol$maxit <- control$maxit # 151217 JN
        mcontrol$trace <- control$trace
# 170924 no longer needed
##        if (control$trace < 1) {mcontrol$iprint <- -1} else {mcontrol$iprint <- control$trace} 
        if (control$trace > 0) cat("lbfgsb3c:control$have.bounds =",control$have.bounds,"\n")
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (control$have.bounds) { 
            slower <- lower/pscale
            supper <- upper/pscale
            ans <- try(lbfgsb3c::lbfgsb3c(par=spar, fn=efn, gr=egr, lower = slower, 
                upper = supper, control=mcontrol)) # explicit pkg in call 170919
            if (! is.numeric(ans$convergence) ) {
                ans$convergence <- 7777 # ?? fixup later
                if (control$trace > 0) {
                   cat("lbfgsb3c - code 7777, non-numeric value of ans$convergence\n")
                   print(ans)
                }
             }
        } else {
            ans <- try(lbfgsb3c::lbfgsb3c(par=spar, fn=efn, gr=egr, control=mcontrol))
        }
        if (! inherits(ans, "try-error")) {
 ## Need to check these carefully. Changed 20191202 for lbfgsb3c !!?
            ans$par <- ans$par*pscale
            ans$niter <- NULL # loss of information
         } else {
            if (control$trace>0) cat("lbfgsb3c failed for current problem \n")
            ans<-list(fevals=NA) # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par<-rep(NA,npar)
            ans$convergence<-9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[1] <- NA
         }
         ## return(ans)
      }  ## end if using lbfgsb3c
## --------------------------------------------
      else if (method == "lbfgs") {# Use unconstrained method from lbfgs package
        if (control$trace > 1) cat("lbfgs\n")
        # Following seems to be needed to avoid unwanted output
        if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
        if (control$trace > 1) cat("lbfgs:control$have.bounds =",control$have.bounds,"\n")
        ans <- list() # to define the answer object
        errmsg <- NA
        class(ans)[1] <- "undefined" # initial setting
        if (control$have.bounds) {
              cat("control$have.bounds seems TRUE\n")
              if (control$trace > 0) cat("lbfgs::lbfgs cannot handle bounds\n")
              errmsg <- "lbfgs::lbfgs cannot handle bounds\n"
            ##  stop("lbfgs::lbfgs tried with bounds")
            class(ans)[1] <- "try-error"            
        }
        if (is.null(egr)) {
            if (control$trace > 0) cat("lbfgs::lbfgs MUST have gradient provided\n")
            errmsg <- "lbfgs::lbfgs MUST have gradient provided\n"
            class(ans)[1] <- "try-error"            
        }
        ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
        if (inherits(ans, "undefined")){ #?? need to explain why no dots??
            ans <- try(lbfgs::lbfgs(efn, egr, vars=spar, 
                    invisible=invisible))
        }
        if (! inherits(ans, "try-error")) {
        ## Need to check these carefully!!?
            ans$par <- ans$par*pscale
            ans$counts[1] <- NA # lbfgs seems to have no output like this
            ans$counts[2] <- NA
         } else {
            if (control$trace>0) cat("lbfgs failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
         ## return(ans)
      }  ## end if using lbfgs
  ## --------------------------------------------
  else if (method == "subplex") {# Use unconstrained method from subplex package
    if (control$trace > 1) cat("subplex\n")
    if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
    if (control$trace > 1) cat("subplex:control$have.bounds =",control$have.bounds,"\n")
    ans <- list() # to define the answer object
    if (control$trace > 0) warning("subplex has no trace mechanism")
    class(ans)[1] <- "undefined" # initial setting
    if (control$have.bounds) {
      cat("control$have.bounds seems TRUE\n")
      if (control$trace > 0) cat("subplex::subplex cannot handle bounds\n")
      stop("subplex::subplex cannot handle bounds")
    }
    ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
    if (class(ans)[1] == "undefined"){
       ans <- try(subplex::subplex(par=spar, fn=efn, control=list(maxit=control$maxfeval)))
#       ans <- try(subplex::subplex(par=spar, fn=efn, control=list(maxit=control$maxfeval), ...))
    }
       if (!inherits(ans, "try-error") && (ans$convergence != -2)) {
       ## Need to check these carefully!!?
       ans$par <- ans$par*pscale
       ans$counts[1] <- ans$count
       ans$counts[2] <- NA
       ans$count <- NULL
       ccode <- ans$convergence
       ans$convergence <- 9999
       ans$message <- paste("subplex:",ans$message)
       if ((ccode == 0) || (ccode == 1)) {
                ans$convergence <- 0 
                if (ccode == 0) { ans$message <- "subplex: success" }
                else { ans$message <- "subplex: Limit of precision reached" }
           } # converged OK
           else {if (ccode == -1) {
                    ans$convergence <- 1
                    ans$message <- "subplex: function evaluation limit reached"
                } # effort limit
       }
    } else { 
      if (ccode == -2) {
            ans$convergence <- 20
      }
      else { ans$convergence <- 9999}
      ans$value <- control$badval
      ans$par <- rep(NA,npar)
      ans$counts[1] <- NA
      ans$counts[2] <- NA
      ans$hessian <- NULL
    }
  }  ## end if using subplex
## --------------------------------------------
  else if (method == "mla") {# Use unconstrained method from marqLevAlg
      if (control$trace > 1) cat("mla\n")
      # Following seems to be needed to avoid unwanted output
      #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
      tans <- list() # to define the answer object
      errmsg <- NA
      class(tans)[1] <- "undefined" # initial setting
      if (control$have.bounds) {
         cat("control$have.bounds seems TRUE\n")
         if (control$trace > 0) cat("marqLevAlg::mla cannot handle bounds\n")
         errmsg <- "marqLevAlg::mla cannot handle bounds\n"
         class(tans)[1] <- "try-error"            
      }
      ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
      if (inherits(tans, "undefined")){
         if (control$trace > 0) pinfo <- TRUE else pinfo <- FALSE
         tans <- try(marqLevAlg::mla(b=spar, fn=efn, gr=egr, hess=ehess,
                     print.info=pinfo))
      }
      if (control$trace > 3) {
        cat("interim answer:")
        str(tans)
      }
      ans <- list() # to ensure set as list
      if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
          ans$par <- tans$b*pscale
          ans$value <- tans$fn.value
          ans$counts[1] <- NA 
          ans$counts[2] <- tans$ni
          ans$convergence <- tans$istop-1
          if (ans$convergence > 2) ans$convergence <- 9999
          tans <- NULL # cleanup
         } else {
            if (control$trace>0) cat("mla failed for current problem \n")
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
      }  ## end if using mla
## --------------------------------------------
  else if (method == "slsqp") {# slsqp from nloptr
      if (control$trace > 1) cat("slsqp\n")
      # Following seems to be needed to avoid unwanted output
      #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
      ans <- list() # to define the answer object
      errmsg <- NA
#      nlalg <- "NLOPT_LD_SLSQP"
      class(ans)[1] <- "undefined" # initial setting
      if (inherits(ans, "undefined")){
         dotstuff<-list(...)
         if (is.null(dotstuff)) { cat("No ...\n") }
          maxeval<-control$maxfeval
          ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
          if (length(slower) == 1) slower<-rep(slower, npar)
          if (length(supper) == 1) supper<-rep(supper, npar)
          tans <- try(nloptr::slsqp(x0=spar, efn, egr, lower=slower, upper=supper,
                        control=list(maxeval=maxeval)) )
      }
      if (control$trace > 3) {
        cat("interim answer:")
        str(tans)
      }
      if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
          ans$par <- tans$par*pscale
          ans$value <- tans$value
          ans$counts[1] <- NA 
          ans$counts[2] <- tans$iter
          if (tans$convergence >= 1) { ans$convergence <- 0 }
          # ?? >1 is supposed to be OK, <0 an error, but 1 not defined.
          else { if (tans$convergence == 0) {ans$convergence <- 9991 }
                 else { ans$convergence <- 9999+tans$convergence }
          if (ans$convergence > 2) { ans$convergence <- 9999 }
          }
          if (ans$convergence > 2) ans$convergence <- 9999
          tans <- NULL # cleanup
         } else {
            if (control$trace>0) cat("slsqp failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
      }  ## end if using slsqp
## --------------------------------------------
  else if (method == "tnewt") {# tnewton from nloptr
      if (control$trace > 1) cat("tnewt\n")
      # Following seems to be needed to avoid unwanted output
      #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
      ans <- list() # to define the answer object
      errmsg <- NA
#      nlalg <- "NLOPT_LD_SLSQP"
      class(ans)[1] <- "undefined" # initial setting
      if (inherits(ans, "undefined")){
         dotstuff<-list(...)
         if (is.null(dotstuff)) { cat("No ...\n") }
          maxeval<-control$maxfeval
          if (length(slower) == 1) slower<-rep(slower, npar)
          if (length(supper) == 1) supper<-rep(supper, npar)
          tans <- try(nloptr::tnewton(x0=spar, efn, egr, lower=slower, upper=supper,
                        control=list(maxeval=maxeval)) )
      }
      if (control$trace > 3) {
        cat("interim answer:")
        str(tans)
      }
      if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
          ans$par <- tans$par*pscale
          ans$value <- tans$value
          ans$counts[1] <- NA 
          ans$counts[2] <- tans$iter
          if (tans$convergence >= 1) { ans$convergence <- 0 }
          # ?? >1 is supposed to be OK, <0 an error, but 1 not defined.
          else { if (tans$convergence == 0) {ans$convergence <- 9991 }
                 else { ans$convergence <- 9999+tans$convergence }
          if (ans$convergence > 2) { ans$convergence <- 9999 }
          }
          if (ans$convergence > 2) ans$convergence <- 9999
          tans <- NULL # cleanup
         } else {
            if (control$trace>0) cat("tnewt failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
      }  ## end if using tnewton
## --------------------------------------------
  else if (method == "anms") {# Use nelder-mead method from pracma
      if (control$trace > 1) cat("anms\n")
      # Following seems to be needed to avoid unwanted output
      #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
      ans <- list() # to define the answer object
      errmsg <- NA
      class(ans)[1] <- "undefined" # initial setting
      ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
      if (inherits(ans, "undefined")){
        if (control$have.bounds) {
              if (control$trace > 0) cat("anms cannot handle bounds\n")
              errmsg <- "anms cannot handle bounds\n"
              stop(errmsg)
              ans <- list()
              class(ans)[1] <- "try-error"
          } else {
              tans <- try(pracma::anms(fn=efn, x0=spar, maxfeval=control$maxfeval))
#              tans <- try(pracma::anms(fn=efn, x0=spar, maxfeval=control$maxfeval, ...))
          }
      }
      if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
          ans$par <- tans$xmin*pscale
          ans$value <- tans$fmin
          ans$counts[1] <- tans$nfeval
          ans$counts[2] <- NA
          ans$convergence<-0
          ans$hessian <- NULL
          if (tans$nfeval >= control$maxfeval) { ans$convergence <- 1 }
          tans <- NULL # cleanup
         } else {
            if (control$trace>0) cat("anms failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
      }  ## end if using anms
## --------------------------------------------
  else if (method == "pracmanm") {# Use nelder_mead from pracma, Gao-Han adaptive NelderMead 
    if (control$trace > 1) cat("pracmanm\n")
    # Following seems to be needed to avoid unwanted output
    #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
    ans <- list() # to define the answer object
    errmsg <- NA
    class(ans)[1] <- "undefined" # initial setting
    ## if(nsctrl > 0) { stop("There are no extra controls set up for ",method) }
    if (inherits(ans, "undefined")){
      if (control$have.bounds) {
        if (control$trace > 0) cat("pracmanm cannot handle bounds\n")
        errmsg <- "pracmanm cannot handle bounds\n"
        stop(errmsg)
        ans <- list()
        class(ans)[1] <- "try-error"
      } else {
        pnmtol <- 1.0e-08 # default in pracma
        if (! is.null(mcontrol$pracmanmtol)) pnmtol <- mcontrol$pracmanmtol
        tans <- try(pracma::nelder_mead(fn=efn, x0=spar, tol=pnmtol, maxfeval=control$maxfeval))
      }
    }
    if (control$trace > 3) {
        cat("interim answer:")
        str(tans)
    }
    if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
      ans$par <- tans$xmin*pscale
      ans$value <- tans$fmin
      ans$counts[1] <- tans$count
      ans$counts[2] <- NA
      ans$convergence<-0
      attr(ans$convergence, "restarts") <- tans$info$restarts
      ans$hessian <- NULL
      ans$message <- tans$errmess
      if (tans$count >= control$maxfeval) { ans$convergence <- 1 }
      tans <- NULL # cleanup
    } else {
      if (control$trace>0) cat("pracmanm failed for current problem \n")
      ans<-list() # ans not yet defined, so set as list
      ans$value <- control$badval
      ans$par <- rep(NA,npar)
      ans$convergence <- 9999 # failed in run
      ans$counts[1] <- NA
      ans$counts[2] <- NA # was [1] until 20211122
      ans$hessian <- NULL
      if (! is.na(errmsg)) ans$message <- errmsg
    }
  }  ## end if using pracmanm
  ## --------------------------------------------
  else if (method == "nlnm") {# neldermead from nloptr
      if (control$trace > 1) cat("nlnm\n")
      # Following seems to be needed to avoid unwanted output
      #  if (control$trace < 1) {invisible <- 1} else {invisible <- 0}
      ans <- list() # to define the answer object
      errmsg <- NA
      class(ans)[1] <- "undefined" # initial setting
      if (inherits(ans, "undefined")){
         dotstuff<-list(...)
         if (is.null(dotstuff)) { cat("No ...\n") }
          maxeval<-control$maxfeval
          if (length(slower) == 1) slower<-rep(slower, npar)
          if (length(supper) == 1) supper<-rep(supper, npar)
          tans <- try(nloptr::neldermead(x0=spar, efn, lower=slower, upper=supper,
                        control=list(maxeval=maxeval)) )
      }
      if (control$trace > 3) {
        cat("interim answer:")
        str(tans)
      }
      if (! inherits(tans, "try-error")) { ## Need to check these carefully!!?
          ans$par <- tans$par*pscale
          ans$value <- tans$value
          ans$counts[1] <- NA 
          ans$counts[2] <- tans$iter
          if (tans$convergence >= 1) { ans$convergence <- 0 }
          # ?? >1 is supposed to be OK, <0 an error, but 1 not defined.
          else { if (tans$convergence == 0) {ans$convergence <- 9991 }
                 else { ans$convergence <- 9999+tans$convergence }
          if (ans$convergence > 2) { ans$convergence <- 9999 }
          }
          if (ans$convergence > 2) ans$convergence <- 9999
          tans <- NULL # cleanup
         } else {
            if (control$trace>0) cat("nlnm failed for current problem \n")
            ans<-list() # ans not yet defined, so set as list
            ans$value <- control$badval
            ans$par <- rep(NA,npar)
            ans$convergence <- 9999 # failed in run
            if (is.null(egr)) ans$convergence <- 9998 # no gradient
            ans$counts[1] <- NA
            ans$counts[2] <- NA # was [1] until 20211122
            ans$hessian <- NULL
            if (! is.na(errmsg)) ans$message <- errmsg
         }
      }  ## end if using nlnm (nloptr::neldermead)
## --------------------------------------------
  ## END OF optimrx extra methods
# ---  UNDEFINED METHOD ---
      else { errmsg<-paste("UNDEFINED METHOD:", method, sep='')
             stop(errmsg, call.=FALSE)
      }
# Exit from routine
      ans$value <- ans$value * fnscale # reset for maximum
      attr(ans$value, "fname") <- fname # add the name if available 20220222
      attr(ans$value, "method") <- method # add method
      attr(ans$par, "status") <- rep("?", npar)
      ans$scounts<-c(optsp$kfn, optsp$kgr, optsp$khe) # counts from function calls
      if ((ans$convergence < 3) && control$have.bounds) { # add indicators !!! Note poss non conv.
         attr(ans$value, "ptype") <- "C"
#         solbds<-bmchk(ans$par, lower = lower, upper = upper, bdmsk = msk,
         solbds<-bmchk(ans$par, lower = lower, upper = upper, bdmsk = NULL,
             trace = control$trace, shift2bound = FALSE)
         attr(ans$par, "status") <- solbds$bchar
      }
      else {
         attr(ans$par, "status") <- rep(" ",npar)
         attr(ans$value, "ptype") <- "U"
      }         
      if (savehess) { # compute hessian
         if (is.null(orig.hess)){
           if (is.null(orig.gr) || is.character(orig.gr) ) {
              hes <- hesf(orig.fn, ans$par, ...) # from numDeriv or pracma??
           } else { 
             hes <- jacf(orig.gr, ans$par, ...) # use Jacobian of gradient
             # 20230613: this approach may give asymmetry of matrix
             hes <- 0.5*(hes + t(hes)) # symmetrize
           }
         } 
         else { # we have a hessian function
           hes <- orig.hess(ans$par, ...) # actual fn
         }
      } else { hes <- NULL } # to ensure it is defined
      ans$hessian <- hes
      ans # last statement of routine
} ## end of optimr
