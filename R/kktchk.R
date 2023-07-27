kktchk <- function(par, fn, gr, hess=NULL, upper=NULL, lower=NULL, maximize=FALSE, control=list(), ...) {
# Provide a check on Kuhn-Karush-Tucker conditions based on quantities
# already computed. Some of these used only for reporting.
##
# Input:
#  par = a single vector of starting values
#  fval = objective function value
#  ngr = gradient evaluated at parameters par
#  nHes = Hessian matrix evaluated at the parameters par
#  nbm = number of active bounds and masks from gHgenb (gHgen returns 0)
#  maximize = logical TRUE if we want to maximize the function. Default FALSE.
#  control = list of controls, currently, 
#            kkttol=1e-3, kkt2tol=1e-6, ktrace=FALSE
#  ... = dot arguments
#
# Output: A list of four elements, namely,
#  gmax = max abs gradient element
#  evratio = ratio of smallest to largest eigenvalue of estimated Hessian
#  kkt1 = logical flag: TRUE if gradient KKT test is satisfied to tolerances
#         otherwise FALSE. Note that decisions are sensitive to scaling and tolerances
#  kkt2 = logical flag: TRUE if Hessian KKT test is satisfied to tolerances
#         otherwise FALSE. Note that decisions are sensitive to scaling and tolerances
# 
# NOTE: Does not do projections for bounds and masks! 
#################################################################
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
#   print(control)
   if (is.null(control$trace)) control$trace <- 0
   if (is.null(control$kkttol)) kkttol<-1e-3 else kkttol<-control$kkttol
   if (is.null(control$kkt2tol)) kkt2tol<-1e-6 else kkt2tol<-control$kkt2tol
   if (control$trace > 0) { cat("kktchk: kkttol=",kkttol,"   kkt2tol=",kkt2tol,
            "  control$trace=",control$trace,"\n") }
   dotargs <- list(...)
#   cat("dotargs:\n")
#   print(dotargs)
   fval <- fn(par, ...) # Note: This requires extra effort to the optimization. 
   # It could be avoided by some clever management of the computations, but there
   # is then the risk we do not get the function value that matches the parameters.
   npar<-length(par)
   if (control$trace > 0) { 
      cat("KKT condition testing\n") 
      cat("Number of parameters =",npar," fval =",fval,"\n") 
   }

   bdout <- bmchk(par, lower = lower, upper = upper, bdmsk = NULL, 
                 trace = control$trace, shift2bound = FALSE) 
   nfree <- sum(bdout$bdmsk[which(bdout$bdmsk==1)])
   nbm <- npar - nfree
   if (control$trace > 0) cat("Number of free parameters =",nfree,"\n")

   kkt1<-NA
   kkt2<-NA
   ngr <- NA
   nHes <- NA
   # test gradient
   if (is.null(gr)) stop("kktchk: A gradient function (or approximation method) MUST be supplied")
   if (is.character(gr)) {
      ngr <- do.call(gr, list(par, fn, ...))
   } else {
      ngr <- gr(par, ...)
   }
   if (control$trace > 0) {
      cat("gradient:")
      print(ngr)
   }
   pngr <- ngr
   pngr[which(bdout$bdmsk != 1)] <- 0.0 # set up projected gradient
   if (control$trace > 0) {
      cat("projected gradient:")
      print(pngr)
   }   
   gmax<-max(abs(pngr)) # need not worry about sign for maximizing
   if (control$trace > 0) {
      cat("max abs projected gradient element =",gmax,"  test tol = ",kkttol*(1.0+abs(fval)),"\n")
   }
   kkt1<-(gmax <= kkttol*(1.0+abs(fval)) ) # Is this sensible?
   if (control$trace > 0) {cat("KKT1 result = ",kkt1,"\n") }

   if (is.null(hess)) { # ?? do we want to tag this so we know what type of computation? 
      if (is.character(gr)){
         nHes <- hessian(func=fn, par, ...) # use numDeriv
      } else {
         nHes <- jacobian(gr, par, ...)
      } 
   } else { 
      nHes <- hess(par, ...)
   }
   if (maximize) {
      nHes<- -nHes
      if (control$trace > 0) cat("Maximizing: use negative Hessian\n")
   }

   # Could provide both free parameter and constrained parameter gradient and Hessian measures
   # Decided to only provide "constrained" ones
#   hev<- try(eigen(nHes)$values, silent=TRUE) # 091215 use try in case of trouble, 
#                                              # 20100711 silent
#   if (!inherits(hev, "try-error")) {
#     if (control$trace > 0) {
#        cat("Hessian eigenvalues of unconstrained Hessian:\n")
#        print(hev) # no check for errors
#     }
#   } else { 
#      warning("Error during computation of unconstrained Hessian eigenvalues") 
#      if(control$trace > 0) cat("Hessian eigenvalue calculation (unconstrained) has failed!\n") 
#       # JN 111207 added
#      hev <- rep(NA, npar) # try to avoid stopping the run, but there is no useful info
#   }

   pHes <- nHes # projected Hessian
#   cat("max asymmetry:", max(abs(pHes-t(pHes))), "\n")
   if (! isSymmetric(unname(pHes))  ) {
      warning("kktchk: pHes not symmetric -- symmetrizing")
      pHes <- 0.5*(pHes + t(pHes))
   }
   pHes[which(bdout$bdmsk != 1), ] <- 0.0
   pHes[ ,which(bdout$bdmsk != 1)] <- 0.0
   if (nfree > 0) {
     phev<- try(eigen(pHes)$values, silent=TRUE) # 091215 use try in case of trouble, 
                                              # 20100711 silent
     if (! inherits(phev,"try-error")) {
        if (control$trace > 0) {
          cat("Hessian eigenvalues of constrained Hessian:\n")
          print(phev) # no check for errors
        }
        # now look at Hessian
        negeig<-(phev[npar] <= (-1)*kkt2tol*(1.0+abs(fval))) # 20100711 kkt2tol
        evratio<-phev[npar-nbm]/phev[1]
        # If non-positive definite, then there will be zero eigenvalues (from the projection)
        # in the place of the "last" eigenvalue and we'll have singularity.
        # WARNING: Could have a weak minimum if semi-definite.
        kkt2<- (evratio > kkt2tol) && (! negeig) 
        if (control$trace > 0) {
          cat("KKT2 result = ",kkt2,"\n") 
        }
        ans<-list(gmax,evratio,kkt1,kkt2, phev, ngatend=ngr, nhatend=nHes)
        names(ans)<-c("gmax","evratio","kkt1","kkt2", "hev", "ngatend", "nhatend")
        return(ans)
     } else {
        warning("Eigenvalue failure for projected Hessian")
        if(control$trace > 0) cat("Hessian eigenvalue calculation (projected) has failed!\n") 
        # JN 111207 added
        phev <- rep(NA, npar) # try to avoid stopping the run, but there is no useful info
        evratio <- NA
        ans<-list(gmax,evratio,kkt1,kkt2, phev, ngatend=ngr, nhatend=nHes)
        names(ans)<-c("gmax","evratio","kkt1","kkt2", "hev", "ngatend", "nhatend")
        return(ans)
     }
  } else {
     warning("All parameters are constrained")
     ans <- list(gmax=0.0, evratio = NA, kkt1=TRUE, kkt2=TRUE, hev=rep(0,npar), ngatend=ngr, nhatend=nHes)
     # this is the fully constrained case
  } # end kkt test
} ## end of kktcchek
