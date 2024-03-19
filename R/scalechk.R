################### scalechk #######################
scalechk<-function(par, lower=lower, upper=upper, bdmsk=NULL, dowarn=TRUE){
   # a function to check the initial parameters and bounds for inputs to optimization codes
   # Arguments:
   #   par -- starting parameters supplied 
   #    lower, upper -- lower and upper bounds supplied
   #    bdmsk -- bounds and masks indicator
   #
   # Returns:
   #   list(lpratio, lbratio) -- the log of the ratio of largest to smallest parameters
   #      and bounds intervals (upper-lower) in absolute value (ignoring Inf, NULL, NA)
   #######################################
   if (is.null(par)) { stop("Null parameter vector") }
   npar<-length(par)
   if (is.null(lower)) { 
      if (dowarn) warning("Null lower bounds vector")
      lower<-rep(-Inf,npar)
   }
   if (is.null(upper)) { 
      if (dowarn) warning("Null upper bounds vector")
      upper<-rep(Inf,npar)
   }
   if (length(lower) == 1) lower<-rep(lower,npar)
   if (length(upper) == 1) lower<-rep(lower,npar)
   # Ignore masked parameters for scaling measures
   if (is.null(bdmsk)) bdmsk <- rep(1,npar) # ensure bdmsk defined
   # Only compute scaling measures for non-masked parameters.
   lower<-lower[which(bdmsk==1)]
   upper<-upper[which(bdmsk==1)]
   newpar<-par[which(bdmsk==1)]
   if ( any(is.infinite(newpar)) ) stop("scalechk:  inadmissible infinite parameter value.")
   else newpar<-abs(newpar) 
   logpar<-log10(newpar[which(newpar > 0)]) # Change 20100711
   # Now check the bounds. Only consider finite ones.
   newlower<-abs(lower[which(is.finite(lower))])
   loglower<-log10(newlower[which(newlower>0)]) # Change 20100711
   newupper<-abs(upper[which(is.finite(upper))])
   logupper<-log10(newupper[which(newupper>0)]) # Change 20100711
   # compute the scalling measures
   bddiff<-upper-lower
   bddiff<-bddiff[which(is.finite(bddiff))] # Only consider finite bounds for scale.
   lbd<-log10(bddiff[which(bddiff > 0)]) # Change 20100711
   if(length(logpar) > 0) ## Thanks to J Laake 140904
     lpratio <- max(logpar) - min(logpar)
   else
     lpratio <- 0 
   if (length(lbd) > 0) {
      lbratio<-max(lbd)-min(lbd)
   } else { 
      lbratio <- NA
   }
   ratios<-list(lpratio=lpratio,lbratio=lbratio) # return(ratios)
}
################### end scalechk #######################
