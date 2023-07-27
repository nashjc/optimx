axsearch<-function(par, fn=NULL, fmin=NULL, lower=NULL, upper=NULL, bdmsk=NULL, 
              control=list(), ...){
# J C Nash 2011-8-2
# Axial search around supposed MINIMUM
# For maximization, user MUST create mfn <- (-1)*fn first.
## Get controls
  npar<-length(par) # number of parameters
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
  ctrl<-NULL # and remove ctrl to avoid trouble
## Set controls we need
  bigval <- control$bigval
  trace <- control$trace
  reltest <- control$reltest
  epst <- control$reltest*control$grtesttol
## BEGIN
  fback<-rep(NA,npar) # for backward step function values
  ffwd <-fback        # for forward step function values
  parstep<- fback     # for the step sizes
  tilt   <- fback     # for the tilts
  roc    <- fback     # for radii of curvature
  bestfn <- Inf  # to record best function found OTHER than fmin (if lower, then restart?)
  par0<-par
# check bounds and masks
  if (is.null(lower)) {
     lower<-rep(-Inf,npar)
  } else {
     if (length(lower) == 1) lower<-rep(lower,npar)
     else if (length(lower) != npar) stop("lower bounds vector of wrong length")
  }
  if (is.null(upper)) {
     upper<-rep(+Inf,npar)
  } else {
     if (length(upper) == 1) upper<-rep(upper,npar)
     else if (length(upper) != npar) stop("upper bounds vector of wrong length")
  }
  if (is.null(bdmsk)) {
     bdmsk<-rep(1,npar)
  } else {
     if (length(bdmsk != npar)) stop("bdmsk bounds and masks vector of wrong length")
  }
## Set fmin
  if (is.null(fmin)) {
     fmin<-fn(par,...) # evaluate the function if we don't have a value
  }
##  fmin0<-rep(fmin,npar) # Needed only for output at the end

  for (j in 1:npar) { # loop over parameters
     if (bdmsk[j] != 0) { # check if NOT masked
  # Parameter may be in bounds, but not ON bound. Is function defined ON bound?
      pkeep<-par[j]
      pstep<-epst*(abs(pkeep+epst))
      parstep[j]<-pstep
      # step backwards
      parj<-pkeep-pstep
      # cat(j, parj, lower[j],"  pstep=",pstep," pkeep=",pkeep,"\n")  
      if (parj < lower[j]) { # out of bounds
         fb<-bigval # set to provide tilt and possibly roc, but leave NA in fback[] 
      } else { if ((reltest+parj)==(reltest+pkeep)) { # no change in parameter
                  fb<-fmin # no change in function
               } else {
                  par[j]<-parj
                  fb<-fn(par,...) # check if admissible -- will be set to bigval otherwise
               }
               fback[j]<-fb
      }
      if (fb < bestfn) bestfn <- fb
      # cat("bestfn=",bestfn,"  fb=",fb,"\n") 
      if (fb < fmin) { # lower value found
          if (trace>0) cat(fb, " *** LOWER ***\n")
          break # to end cycle -- parameters are reset at moment
      }
      parj<-pkeep+pstep # step forward
      if (parj > upper[j]) { # out of bounds
         ff<-bigval # set to provide tilt and possibly roc, but leave NA in ffwd[] 
      } else { if ((reltest+parj)==(reltest+pkeep)) { # no change in parameter
                  ff<-fmin # no change in function
               } else {
                  par[j]<-parj
                  ff<-fn(par,...) # check if admissible -- will be set to bigval otherwise
               }
               ffwd[j]<-ff
      }
      if (ff < bestfn) bestfn <- ff
      if (ff < fmin) { # lower value found
         if (trace>0) cat(ff, " *** LOWER ***\n")
         break # to end cycle -- parameters are reset at moment
      }
      par[j]<-pkeep # reset parameter value
      # compute tilts and curvature
      c1<-0.5*(ff-fb)/pstep # linear term -- should be zero
      c2<-(ff+fb-2*fmin)/(2*pstep*pstep) # quadratic term
      if (trace>0) cat("B[",j,"]=",pkeep," pstep=",pstep,"  f-, f+:",fb, ff,"\n")
      c0 <- 1+c1*c1 # denominator for curvature
      c2=c2/(c0*sqrt(c0))
      c0<-bigval # set large in case of singularity
      if (c1 !=0 ) c0 <- 1/c2 #  radius of curvature
      roc[j]<-c0
      tilt[j] <- -45*atan(c1)/atan(1) # REM tilt 
   } # end if masked
     else { # action when masked
        tilt[j]<-0
        roc[j]<-0
        fback[j]<- bigval 
        ffwd[j] <- bigval
  }
} # end of loop over parameters
if (bestfn >= fmin) { bestfn<-fmin } # set it to be sure when we have no better fn
result<-list(bestfn=bestfn, par=par, details=data.frame(par0, fback, rep(fmin, npar), ffwd, parstep, tilt, roc))
} # end axsearch.R

