cnvtst <- function  (alpha, pnorm, xnorm, 
		 dif, ftest, gnorm, gtp, f, flast, g, ipivot, accrcy) {
## Note: trace is integer (changed 180329)
## ---------------------------------------------------------
##  test for convergence
## ---------------------------------------------------------
##  set up
## ---------------------------------------------------------
conv   <- 0;
eps <- .Machine$double.eps
toleps <- sqrt(accrcy) + sqrt(eps);
rtleps <- accrcy + eps;
imax   <- 0;
ltest  <- (flast - f <= -0.5*gtp);
## ---------------------------------------------------------
##  if anti-zigzag test satisfied, test multipliers;
##  if appropriate, modify active set
## ---------------------------------------------------------
   if ( ! ltest) {
      ind   <- which( (ipivot != 0)  & (ipivot != 2))
      if ( length(ind) > 0 ) {  ## how to ensure ind not empty
         t <- -sum(ipivot[ind]*g[ind])
         cmax <- min(t) # why [cmax, imax] = min(t) ?
         imax <- cmax
         if (cmax >= 0) { imax <- 0 }
      }
   }
   if (imax != 0) {
      ipivot[ ind[imax] ] <- 0;
      flast <- f;
   } else {
      conv = ( ( (alpha*pnorm < toleps*(1 + xnorm) )
           && (abs(dif) < rtleps*ftest)
           && (gnorm < accrcy^(1/3)*ftest) ) || 
            ( gnorm < .01*sqrt(accrcy)*ftest) )
   }
   result <- list(conv=conv, flast1=flast, ipivot1=ipivot)
} 

crash <- function (x, low, up) {
##---------------------------------------------------------
## this initializes the constraint information, and
## ensures that the initial point satisfies 
##      low <= x <= up.
## the constraints are checked for consistency.
##---------------------------------------------------------
ierror <- 0 
if (any(low > up)) { ierror = - max(which(low > up))  } 
   # above is check on error in bounds specification
   xnew  <- pmax (low, x) # force params into bounds
   xnew  <- pmin (up, xnew) # No diagnostic!!?
   # we output revised parameters
   n <- length(x)
   ind <- which(low == x) 
   ipivot <- rep(0, n) ## zeros(size(x)) in Matlab
   if (length(ind) > 0) { ipivot[ind] <- -1 }
   ind  <- which(x == up) 
   if (length(ind) > 0) { ipivot[ind] <-  1 }
   ind <- which(low == up) 
   if (length(ind) > 0) { ipivot[ind] <- 2 }
   list(ipivot=ipivot, ierror=ierror, xnew=xnew)
}


gtims <- function (v, x, g, accrcy, xnorm, sfun, ...) {
##---------------------------------------------------------
## compute the product of the Hessian times the vector v;
## store result in the vector gv 
## (finite-difference version)
##---------------------------------------------------------
##   cat("v, x, g:")
##   print(v)
##   print(x)
##   print(g)
   delta <- sqrt(accrcy)*(1 + xnorm)/sqrt(sum(v^2))
   hg <- x + delta*v
   tresult<-sfun(hg, ...)
##   cat("tresult: ")
##   print(tresult)
   gv <- attr(tresult, "gradient")
   gv <- (gv - g)/delta
   gv 
}

initpc <- function(d, upd1, ireset) {
##---------------------------------------------------------
## initialize the diagonal preconditioner (d -- a vector!)
## ---------------------------------------------------------
## global hyk sk yk sr yr yksk yrsr # do we have these?
#  global vectors hyk sk yk sr yr & scalars yksk yrsr
## ---------------------------------------------------------
#%   cat("In initpc -- sk and hyk if upd1 false: ", upd1,"\n")
#%   print(envjn$sk)
   if (upd1) { 
      td <- d
   } else {
      if (ireset) {
         envjn$hyk  <- d * envjn$sk # vector * vector by element?
         sds  <- as.numeric(crossprod(envjn$sk, envjn$hyk)) # matrix multiply
         if (all(envjn$hyk == 0) && trace > 1) { cat("INITPC: envjn$hyk = 0 \n") }
         if (sds == 0 && trace > 1) { cat("INITPC: sds = 0 \n") }
         td   <- d - d*d*envjn$sk*envjn$sk/sds + envjn$yk*envjn$yk/envjn$yksk 
         # by element 
      } else {
         envjn$hyk  <- d * envjn$sr # vector * vector
         sds  <- as.numeric(crossprod(envjn$sr, envjn$hyk))
         srds <- as.numeric(crossprod(envjn$sk, envjn$hyk))
         yrsk <- as.numeric(crossprod(envjn$yr, envjn$sk))
         envjn$hyk  <- d*envjn$sk - envjn$hyk*srds/sds + envjn$yr*yrsk/envjn$yrsr
         td   <- d - d*d*envjn$sr*envjn$sr/sds+envjn$yr*envjn$yr/envjn$yrsr
         sds  <- as.numeric(crossprod(envjn$sk, envjn$hyk))
         td   <- td - envjn$hyk*envjn$hyk/sds + envjn$yk*envjn$yk/envjn$yksk
      }
   }
#%   cat("td:")
#%   print(td)
   ans<-list(td=td)
}


lin1 <- function(p, x, f, alpha, g, sfun, ...){
   ## ---------------------------------------------------------
   ##  line search (naive)
   ## ---------------------------------------------------------
   ##  set up
   ## ---------------------------------------------------------
#   cat("lin1: alpha=",alpha,"  p:\n")
#   print(p)

   if (is.null(alpha)) alpha <- 0

   ierror <- 3 
   xnew   <- x 
   fnew   <- f 
   gnew   <- g 
   maxit  <- 15 
   if (alpha == 0) {
      ierror <- 0
      maxit <- 1
   }
   alpha1 <- alpha 
   ## ---------------------------------------------------------
   ##  line search
   ## ---------------------------------------------------------
   for (itcnt in 1:maxit) {
      xt <- x + alpha1*p 
      fg <- sfun(xt, ...) # Note: added dots 140902
      ft<-fg
      gt<- attr(fg,"gradient") # may simplify later
      if (ft < f) {
         ierror <- 0
         xnew   <- xt
         fnew   <- ft
         gnew   <- gt
         ## cat("about to break in lin1\n")
         break
      }
      alpha1 <- alpha1 / 2
   }
   if (ierror == 3) { alpha1 <- 0 }
   nf1 <- itcnt 

   ## never used in SGN code
   ## if (nargout == 7) { ## What is nargout?
   ##  dfdp <- as.numeric(crossprod(gt, p) )
   ##  varargout{1} <- dfdp ## What is varargout?
   ## end
   result<-list(xnew=xnew, fnew=fnew, gnew=gnew, nf1=nf1,
          ierror=ierror, alpha1=alpha1)
}

lmqnbc <- function (x, sfun, lower, upper, maxit, maxfun, stepmx, accrcy, trace, ...) {
## ---------------------------------------------------------
##  This is a bounds-constrained truncated-newton method.
##  The truncated-newton method is preconditioned by a 
##  limited-memory quasi-newton method (computed by
##  this routine) with a diagonal scaling (routine ndia3).
##  For further details, see routine tnbc.
## ---------------------------------------------------------
##  global hyk sk yk sr yr yksk yrsr
## ---------------------------------------------------------
##  check that initial x is feasible and that the bounds 
##  are consistent
## ---------------------------------------------------------
   n <- length(x)

# JN: Define "globals" here
   gtn<-list(yrsr=0, yksk=0, yr = rep(0, n), yk = rep(0, n), 
        sr = rep(0, n),  sk = rep(0, n),
        hg=rep(0,n), hyk=rep(0,n), hyr=rep(0,n) )
   envjn<<-list2env(gtn) # Note envjn is in parent environment
# end globals

##   [ipivot, ierror, x] = crash(x, lower, upper) 
   if (trace > 1) cat("lmqnbc -- crout:")
   crout<-crash(x, lower, upper)
   if (trace) print(crout)
   ierror <- crout$ierror
   ipivot <- crout$ipivot
   x<-crout$xnew # in case x changed by bounds
   f <- 0 
   g <- rep(0,n) 
   if (ierror != 0) {
      stop('LMQNBC: terminating (no feasible point)')
      ##   return 
   }
## ---------------------------------------------------------
##  initialize variables, parameters, and constants
## ---------------------------------------------------------
   options(digits=5)
   if (trace > 1) cat("stempmx, accrcy, maxfun:",stepmx, accrcy, maxfun, "\n")
   if (trace > 0) cat('  it     nf     cg           f             |g|\n')
   eps <- .Machine$double.eps
   upd1   <- TRUE
   ncg    <- 0 
   conv   <- FALSE
   xnorm  <- max(abs(x))
   ierror <- 0 

   if ( (stepmx < sqrt(accrcy)) || (maxfun < 1) ) { 
      ierror <- -1 
      xstar <- x  
      almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                 nfngr=ncg)
     if (trace > 1) {
        cat("Exiting lmqnbc - almgn:")
        print(almqn)
     }
     return(almqn)  # return 
   }
## ---------------------------------------------------------
##  compute initial function value and related information
## ---------------------------------------------------------
   if (trace > 1) cat("Try initial fn\n")
   fg <- sfun(x, ...)
   nf     <- 1 
   nit    <- 0 
   g<- attr(fg,"gradient")
   f<-fg
   flast  <- f
#   if (is.null(g) ) { ## 160922 change
   if (! is.numeric(g) ) { # try fix 180328
     gnorm <- 1.0/eps 
   } else { gnorm  <- max(abs(g)) } ##  norm(g,'inf') 
## ---------------------------------------------------------
##  Test if Lagrange multipliers are non-negative.
##  Because the constraints are only bounds, the Lagrange
##  multipliers are components of the gradient.
##  Then form the projected gradient.
## ---------------------------------------------------------
   ind <- which((ipivot != 2) & 
             (as.numeric(crossprod(ipivot,g)) >0 ) ) 
   if (length(ind) > 0) {  
      ipivot[ind] <- rep(0, length(ind)) 
   } 
   g <- ztime (g, ipivot) 
#   if (is.null(g) ) { ## 160922 change
   if (! is.numeric(g) ) { # try fix 180328
     gnorm <- 1.0/eps 
   } else { gnorm  <- max(abs(g)) } ##  norm(g,'inf') 
   if (trace > 0) cat(nit,"\t", nf,"\t", ncg,"\t", f,"  ", gnorm,"\n")
   if (trace > 1) {
      print(x)
      print(g)
   }
## ---------------------------------------------------------
##  check if the initial point is a local minimum.
## ---------------------------------------------------------
   ftest <- 1 + abs(f) 
   if (gnorm < .01*sqrt(eps)*ftest) {
      xstar <- x 
      almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
              nfngr=ncg)
      return(almqn)
   }

   if (trace > 1) cat("before set initial, flast=",flast,"\n")

## ---------------------------------------------------------
##  set initial values to other parameters
## ---------------------------------------------------------
   icycle <- n-1 
   ireset <- 0 
   bounds <- TRUE 
   difnew <- 0 
   epsred <- .05 
   fkeep  <- f 
   d      <- rep(1,n) 
## ---------------------------------------------------------
##  ..........main iterative loop..........
## ---------------------------------------------------------
##  compute the search direction
## ---------------------------------------------------------
   argvec <- c(accrcy, gnorm, xnorm) 
##[p, gtp, ncg1, d] <- ...
##	modlnp (d, x, g, maxit, upd1, ireset, bounds, ipivot, ##   argvec, sfun)
   mres  <- modlnp (d, x, g, maxit, upd1, ireset,
             bounds, ipivot, argvec, sfun, ...) 
   ncg1 <- mres$ncg1
   gtp  <- mres$gtp
   d    <- mres$dnew
   p    <- mres$p
   ## cat("p:")
   ## print(p)
   ## tmp<-readline("cont.")

   ncg <- ncg + ncg1 
   while (!conv) {
      oldg <- g 
      pnorm <- max(abs(p)) # norm2(p, 'inf') 
      oldf <- f 
## ---------------------------------------------------------
##  line search
## ---------------------------------------------------------
      pe <- pnorm + eps 
      spe <- stpmax (stepmx, pe, x, p, ipivot, lower, upper) 
      ## cat("spe=",spe,"\n")
      alpha <- step1 (f, gtp, spe) 
      alpha0 <- alpha 
      ## cat("alpha0=",alpha0,"\n")
##   [x_new, f_new, g_new, nf1, ierror, alpha] <- lin1 (p, x,
##       f, alpha0, g, sfun) 
      reslin <- lin1 (p, x, f, alpha0, g, sfun, ...) 
      ierror <- reslin$ierror
      alpha <- reslin$alpha1
## ---------------------------------------------------------
      if ((alpha == 0) && (alpha0 != 0) || (ierror == 3)){ 
          cat('Error in Line Search\n') 
          cat('    ierror = ', ierror, "\n") 
          cat('    alpha  = ',alpha, "\n") 
          cat('    alpha0 = ', alpha0, "\n") 
          cat('    gtp    = ', gtp, "\n") 
          ## ############################
          cat('    |g|     = ', norm2(g), "\n") 
          cat('    |p|     = ', norm2(p), "\n") 
#          tmp <- readline('Hit any key to continue')
          ## ############################
      } 
      ## #######################
      x <- reslin$xnew # need fixup
      f <- reslin$fnew
      g <- reslin$gnew
      nf1 <- reslin$nf1
      nf  <- nf  + nf1 
      nit <- nit +   1 
## ---------------------------------------------------------
##  update active set, if appropriate
## ---------------------------------------------------------
      newcon <- FALSE 
      ## cat("Check active set - alpha, spe, f:",alpha, spe,f,"\n")
      if (abs(alpha-spe) <= 10*eps) {
         newcon <- TRUE
         ierror <- 0 
          ## cat("flast:",flast,"\n")
#         [ipivot, flast] <- modz (x, p, ipivot, lower, upper, 
#             flast, f, alpha) 
         ## cat("ipivot before and after update by modz:\n")
         ## print(ipivot)
         modzres<-modz(x, p, ipivot, lower, upper, 
                      flast, f, alpha) 
         ipivot <- modzres$ipivot1
         flast<-modzres$flast1
         ## print(ipivot)
      }
      if (ierror == 3) { 
         xstar <- x  
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                 nfngr=ncg)
         if (trace > 0) {
            cat("Error updating active set -- almqn:")
            print(almqn)
         }
         return(almqn)
      }
## ---------------------------------------------------------
##  stop if more than maxfun evaluations have been made
## ---------------------------------------------------------
      if (nf > maxfun) { 
         ierror <- 2 
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                     nfngr=ncg)
         if (trace > 0) {
           cat("Too many function evaluations -- almqn:")
           print(almqn)
         }
         return(almqn)
      } 
## ---------------------------------------------------------
##  set up for convergence and resetting tests
## ---------------------------------------------------------
      difold <- difnew 
      difnew <- oldf - f # scalar 
      if (icycle == 1) {
         if (difnew >  2*difold) { epsred <-  2*epsred } 
         if (difnew < .5*difold) { epsred <- .5*epsred } 
      } 
      gv    <- ztime (g, ipivot) 
      gnorm <- max(abs(gv)) 
      ftest <- 1 + abs(f) 
      xnorm <- max(abs(x))
############### DISPLAY ############## 
      if (trace > 0) cat(nit,"\t", nf,"\t", ncg,"\t", f,"  ", gnorm,"\n")
      ## print(x)
      ## print(g)
## ---------------------------------------------------------
##  test for convergence
## ---------------------------------------------------------
##   [conv, flast, ipivot] <- cnvtst (alpha, pnorm, xnorm, ...
##	    difnew, ftest, gnorm, gtp, f, flast, g, ...
##	    ipivot, accrcy) 
      ## cat("before cnvtst, flast=",flast,"\n")
      ctres <- cnvtst (alpha, pnorm, xnorm, 
	       difnew, ftest, gnorm, gtp, f, flast, g, 
	       ipivot, accrcy) 
      conv <- ctres$conv
      flast <- ctres$flast1
      ipivot <- ctres$ipivot1
      ## cat("after cnvtst - conv, flast, ipivot:", conv, flast,"\n")
      ## print(ipivot)
      if (conv) {
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
                     nfngr=ncg)
         return(almqn)
      } 
      g <- ztime (g, ipivot) 
## ---------------------------------------------------------
##  modify data for LMQN preconditioner
## ---------------------------------------------------------
      if (! newcon) { ## 0 is FALSE and 1 is TRUE supposedly
         envjn$yk <- g - oldg 
         envjn$sk <- alpha*p 
         envjn$yksk <- as.numeric(crossprod(envjn$yk, envjn$sk)) 
         ## cat("yksk:",envjn$yksk,"\n")
         ireset <- ( (icycle == n-1) |  
                 (difnew < epsred*(fkeep-f)) ) 
         if (! ireset) { 
            envjn$yrsr <- 
                 as.numeric(crossprod(envjn$yr,envjn$sr)) 
            ireset <- (envjn$yrsr <= 0) 
         } 
         upd1 <- (envjn$yksk <= 0) 
      }
      ## cat("newcon, upd1:", newcon, upd1,"\n")
## ---------------------------------------------------------
##  compute the search direction
## ---------------------------------------------------------
      argvec <- c(accrcy, gnorm, xnorm) 
##   [p, gtp, ncg1, d] <- ...
##	   modlnp (d, x, g, maxit, upd1, ireset, bounds, 
##                   ipivot, argvec, sfun) 
      mres  <- modlnp (d, x, g, maxit, upd1, ireset,
                bounds, ipivot, argvec, sfun, ...) 
      ncg1 <- mres$ncg1
      gtp  <- mres$gtp
      d    <- mres$dnew
      p    <- mres$p

   ## cat("New p:")
   ## print(p)
   ## tmp<-readline("cont.")

      ncg <- ncg + ncg1 
## ---------------------------------------------------------
##  update LMQN preconditioner
## ---------------------------------------------------------
      if (! newcon) { 
         if (ireset) { 
            envjn$sr  <- envjn$sk 
            envjn$yr  <- envjn$yk 
            fkeep  <- f 
            icycle <- 1 
         } else {
            envjn$sr     <- envjn$sr + envjn$sk 
            envjn$yr     <- envjn$yr + envjn$yk 
            icycle <- icycle + 1 
         }
      }
   } # end while !conv
}


lmqn <- function (x, sfun, maxit, maxfun, stepmx, accrcy, trace, ...) {
## ---------------------------------------------------------
##  truncated-newton method for unconstrained minimization
##  (customized version)
## ---------------------------------------------------------
#  global vectors hyk sk yk sr yr & scalars yksk yrsr
## ---------------------------------------------------------
##  set up
## ---------------------------------------------------------
## format compact
## format short e
  n<-length(x)
  # JN: Define globals here. Is it necessary to set up.
   gtn<-list(yrsr=0, yksk=0, yr = rep(0, n), yk = rep(0, n), sr = rep(0, n),  sk = rep(0, n))
   envjn<<-list2env(gtn)
# end globals
   eps <- .Machine$double.eps
   upd1 <- 1 
   ncg  <- 0 
   xnorm  <- max(abs(x)) # norm(x,'inf') 
   ierror <- 0 
   if (stepmx < sqrt(accrcy) || maxfun < 1) {
     ierror <- -1 
     xstar <- x  
     almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror,
              nfngr=ncg)
     return(almqn)
   }
## ---------------------------------------------------------
##  compute initial function value and related information
## ---------------------------------------------------------
   fg <- sfun(x, ...)
#%    print(fg)
   g<- attr(fg, "gradient")
#%   print(g)
   if (is.null(g)) stop("Must have gradient defined for Rtnmin")
   f<-fg
#   if (is.null(g) ) { ## 160922 change
   if (! is.numeric(g) ) { # try fix 180328
     gnorm <- 1.0/eps 
   } else { gnorm  <- max(abs(g)) } ##  norm(g,'inf') 
   nf     <- 1 
   nit    <- 0 
   if (trace)  cat("Itn ",nit," ",nf," ",ncg, " ",f, " ", gnorm,"\n")
## ---------------------------------------------------------
##  check for small gradient at the starting point.
## ---------------------------------------------------------
   ftest <- 1 + abs(f) 
   if (gnorm < .01*sqrt(eps)*ftest) {
     ierror <- 0 
     xstar <- x 
     almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
     return(almqn)
   }
## ---------------------------------------------------------
##  set initial values to other parameters
## ---------------------------------------------------------
   n      <- length(x) 
   icycle <- n-1 
   toleps <- sqrt(accrcy) + sqrt(eps) 
   rtleps <- accrcy + eps 
   difnew <- 0 
   epsred <- .05 
   fkeep  <- f 
   conv   <- FALSE 
   ireset <- 0 
   ipivot <- 0 
   if (trace > 1) cat("end initial values\n")
## ---------------------------------------------------------
##  initialize diagonal preconditioner to the identity
## ---------------------------------------------------------
   d <- rep(1,n) # as a vector 
## ---------------------------------------------------------
##  ..........main iterative loop..........
## ---------------------------------------------------------
##  compute search direction
## ---------------------------------------------------------
   argvec <- c(accrcy, gnorm, xnorm) 
#   cat("call modlnp\n")
   mres  <- modlnp (d, x, g, maxit, upd1, ireset, bounds=FALSE, ipivot, argvec, sfun, ...) 
   p <- mres$p
#   cat("p from first call to modlnp\n")
#   print(p)
#   tmp<-readline("cont.")

   gtp <- mres$gtp
   ncg1<-mres$ncg1
   d <- mres$dnew

#%    cat("d from modlnp:")
#%    print(d)
   ncg <- ncg + ncg1 
   while ( ! conv) { 
#%       cat("top while ")
#%       tmp <- readline("Top of iteration")
      oldg   <- g 
      pnorm  <- max(abs(p)) # norm(p,'inf') 
      oldf   <- f 
## ---------------------------------------------------------
##  line search
## ---------------------------------------------------------
      if (trace > 1) cat("start linesearch\n")
      pe     <- pnorm + eps 
      spe    <- stepmx/pe 
      if (trace > 1) cat("gtp, spe:", gtp, spe,"\n")
      alpha0 <- step1 (f, gtp, spe) 
      reslin <- lin1 (p, x, f, alpha0, g, sfun, ...) 
##      [x, f, g, nf1, ierror, alpha] <- 
      x <- reslin$xnew # need fixup
      f <- reslin$fnew
      g <- reslin$gnew
      nf1 <- reslin$nf1
      ierror <- reslin$ierror
      alpha <- reslin$alpha1
#      cat("after lin1, alpha=",alpha,"\n")
#   tmp<-readline("cont.")

     nf <- nf + nf1 
## ---------------------------------------------------------
      nit <- nit + 1
#   if (is.null(g) ) { ## 160922 change
   if (! is.numeric(g) ) { # try fix 180328
        gnorm <- 1.0/eps 
      } else { gnorm  <- max(abs(g)) } ##  norm(g,'inf') 
### Display info
      if (trace > 0) cat("Itn ",nit," ",nf," ",ncg, " ",f, " ", gnorm,"\n")
      if (ierror == 3) { 
         if (length(ncg) == 0) { ncg <- 0 } # is.null(ncg)?
         xstar <- x 
         almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
         return(almqn)
      }
## ---------------------------------------------------------
##  stop if more than maxfun evalutations have been made
## ---------------------------------------------------------
      if (nf >= maxfun) {
        ierror <- 2 
        xstar <- x 
        almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
        return(almqn)
      }
## ---------------------------------------------------------
##  set up for convergence and resetting tests
## ---------------------------------------------------------
      ftest  <- 1 + abs(f) 
      xnorm  <- max(abs(x)) # norm(x,'inf') 
      difold <- difnew 
      difnew <- oldf - f 
      envjn$yk     <- g - oldg 
      envjn$sk     <- alpha*p 
      if (icycle == 1) {
          if (difnew >   2*difold) { epsred <-   2*epsred }
          if (difnew < 0.5*difold) { epsred <- 0.5*epsred }
      }
## ---------------------------------------------------------
##  convergence test
## ---------------------------------------------------------
      conv <- ( ( (alpha*pnorm < toleps*(1 + xnorm)) &&
                 (abs(difnew) < rtleps*ftest)  &&
                 (gnorm < accrcy^(1/3)*ftest) )||
                 ( gnorm < .01*sqrt(accrcy)*ftest ) )
      if (conv) {
        ierror <- 0 
        xstar <- x 
        almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
        return(almqn)
      }
## ---------------------------------------------------------
##  update lmqn preconditioner
## ---------------------------------------------------------
      envjn$yksk <- as.numeric(crossprod(envjn$yk, envjn$sk) )
      ireset <- ((icycle == n-1) || (difnew < epsred*(fkeep-f)) )
      if ( ! ireset) { 
          envjn$yrsr <- as.numeric(crossprod(envjn$yr, envjn$sr))
          ireset <- (envjn$yrsr <= 0) 
      }
      upd1 <- (envjn$yksk <= 0) 
## ---------------------------------------------------------
##  compute search direction
## ---------------------------------------------------------
      argvec <- c(accrcy, gnorm, xnorm) 
##      cat("ireset =", ireset," upd1=",upd1,"\n")
#%       cat("New d to modlnp:")
#%       print(d)
##   [p, gtp, ncg1, d] <- ...
      mres <- modlnp (d, x, g, maxit, upd1, ireset, 0, ipivot, argvec, sfun, ...) 
      p <- mres$p
      gtp <- mres$gtp
      ncg1<-mres$ncg1
      d <- mres$dnew
      
#%       cat("returned dnew, envjn$john: ", envjn$john)
#%       print(d)
      ncg <- ncg + ncg1 
## ---------------------------------------------------------
##  store information for lmqn preconditioner
## ---------------------------------------------------------
      if (ireset) {
          envjn$sr <- envjn$sk 
          envjn$yr <- envjn$yk 
          fkeep <- f 
          icycle <- 1 
      } else {
          envjn$sr <- envjn$sr + envjn$sk 
          envjn$yr <- envjn$yr + envjn$yk 
          icycle <- icycle + 1 
      }
   } # end while
## [xstar, f, g, ierror] = ..
   almqn<-list(xstar=xstar, f=f, g=g, ierror=ierror, nfngr=ncg)
} # end lmqn


modlnp <- function(d, x, g, maxit, upd1, ireset, bounds, 
         ipivot, argvec, sfun, ...) {
##---------------------------------------------------------
## this routine performs a preconditioned conjugate-gradient
## iteration to solve the Newton equations for a search
## direction for a truncated-newton algorithm. 
## When the value of the quadratic model is sufficiently 
## reduced, the iteration is terminated.
##---------------------------------------------------------
## parameters
##
## p           - computed search direction
## g           - current gradient
## gv,gz1,v    - scratch vectors
## r           - residual
## d           - diagonal preconditoning matrix
## feval       - value of quadratic function
##------------------------------------------------------------
## initialization
##------------------------------------------------------------

   if (is.null(d)) stop("Null d")
   ## print(x)
   ## print(g)
   ## cat(bounds,"\n")

   accrcy <- argvec[1] 
   gnorm  <- argvec[2] 
   xnorm  <- argvec[3] 

   if (maxit == 0) {
      p    <- -g 
      gtp  <- as.numeric(crossprod(p, g))
      ncg1 <- 1
      dnew <- d
      if (sqrt(sum(p^2))==0) {
         ## cat("MODLNP 01: |p| = 0\n") 
         ## pause(1) 
      }
#%       cat("modlnp - dout01:")
#%       print(dnew)
      result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
      return(result) 
   }

   first <- 1
   tol   <- 1e-6   ######### was 1.d-12  #########
   qold  <- 0
   ainit  <- initpc (d, upd1, ireset)
   dnew<-ainit$td
   
#%    cat("after initpc dnew:")
#%    print(dnew)
   r     <- -g
   v     <- rep(0, length(r))
   p     <- v
   gtp   <- as.numeric(crossprod(p, g))

   rho  <- rep(0, maxit+1)
   beta <- rep(0, maxit)
   v.gv <- rep(0, maxit)

   rho[[1]] <- as.numeric(crossprod(r))

##------------------------------------------------------------
## main iteration (conjugate-gradient iterations for Ax = b)
##------------------------------------------------------------
   ind <- 0
   ncg1 <- 0
   for (k in 1:maxit) {
      ncg1 <- ncg1 + 1
      if (bounds) { r  <- ztime(r, ipivot) }
      amsolve <- msolve (r, upd1, ireset, first, d) 
      zk<-amsolve$y
      
      if (bounds) { zk <- ztime (zk, ipivot) }
#%       cat("r:")
#%       print(r)
#%       cat("zk:")
#%       print(zk)
      rz <- as.numeric(crossprod(r,zk))  

      if (rz/gnorm < tol) {
         ind <- 80 
         if (sqrt(sum(p^2))==0) {
            p <- -g
            gtp <- as.numeric(crossprod(p, g))
         }
#%          cat("modlnp - dout - ind80:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      } 
      if (k > 1) {
         beta[k] <- rz/rzold 
      } else {
         beta[k] <- 0 
      }
#%       cat("beta[",k,"] =",beta[k],"\n")
      v <- zk + beta[k]*v 
      if (bounds) { v  <- ztime( v, ipivot) }
      ## cat("about to call gtims\n")
      gv <- gtims(v, x, g, accrcy, xnorm, sfun, ...) 
#%     cat("After gtims: gv=")
#%     print(gv)
      ## cat(bounds,"\n")
      if (bounds) { gv <- ztime (gv, ipivot) }
#%       cat("gv, v:")
      ## cat("gv=", gv, "\n")
#%       print(gv)
#%       print(v)
      v.gv[[k]] <- as.numeric(crossprod(v, gv))  ##  v.gv has underscore!!?
#%       cat("v.gv[[",k,"]]=",v.gv[[k]],"\n")
      if (v.gv[[k]]/gnorm < tol) { 
         ind <- 50  
         if (sqrt(sum(p^2))==0) {
            p <- -g ## Mod SGN 140912
            if (bounds) {
              p <- ztime(p, ipivot)
            }
            gtp <- crossprod(p, g)
##            cat("MODLNP 03: |p| = 0 \n") 
            ## pause(1) 
         }
#%          cat("modlnp - dout03:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      }
      dnew <- ndia3(dnew, v, gv, r, v.gv[[k]]) # NB need to return something but it doesn't get used
#%       cat("after ndia3 below dout03, dnew:")
#%       print(dnew)
##------------------------------------------------------------
## compute current solution and related vectors
##------------------------------------------------------------
      alpha <- rz / v.gv[[k]]
      p <- p + alpha* v 
      r <- r - alpha*gv 

      rho[[k+1]] <- as.numeric(crossprod(r))

##------------------------------------------------------------
## test for convergence
##------------------------------------------------------------
      gtp <- as.numeric(crossprod(p,g))
      pr <-  as.numeric(crossprod(r,p))
      q <- (gtp + pr) / 2 
      qtest <- k * (1 - qold/q) 
      if (qtest <= 0.5)  {
         if (sqrt(sum(p^2))==0) {
##            cat("MODLNP 04: |p| = 0\n") 
            ## pause(1) 
         }
#%          cat("modlnp - dout04:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      }
      qold <- q 
##------------------------------------------------------------
## perform cautionary test
##------------------------------------------------------------
      if (gtp > 0) {
         ind <- 40  
         if (sqrt(sum(p^2))==0) {
##            cat("MODLNP 05: |p| = 0 \n") 
            ##   pause(1) 
         }
#%          cat("modlnp - dout05:")
#%          print(dnew)
         result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
         return(result) 
      } 
      rzold <- rz 
   } ## end loop over k
   k <- k-1 
##------------------------------------------------------------
## terminate algorithm
##------------------------------------------------------------
   if (ind == 40) {
      p   <- p - alpha*v 
   } 
##------------------------------------------------------------
   if (ind == 50 && k <= 1) {
      amsolve <- msolve (g, upd1, ireset, first, d) 
      p<-amsolve$y
      
      p <- -p 
      if (bounds) { p <- ztime (p, ipivot) }
   }
##------------------------------------------------------------
   if (ind == 80 && k <= 1) {
      p <- -g 
      if (bounds) { p <- ztime (p, ipivot) }
   }
##------------------------------------------------------------
## store new diagonal preconditioner
##------------------------------------------------------------
   gtp  <- as.numeric(crossprod(p, g)) 
   ncg1 <- k + 1 
   if (sqrt(sum(p^2))==0) { 
##       cat("MODLNP 06: |p| = 0 \n")
       ##    pause(1) 
   }
   ## cat("modlnp - dout06:")
#%    print(dnew)
   result<-list(p=p, gtp=gtp, ncg1=ncg1, dnew=dnew)
   return(result) 
}

modz <- function (x, p, ipivot, low, up, flast, f, alpha){
##---------------------------------------------------------------------
## update the constraint matrix if a new constraint is encountered
##---------------------------------------------------------------------
   eps <- .Machine$double.eps
   indl <- which(ipivot == 0 & p < 0)
   if (length(indl) > 0) {
      toll <- 10 * eps * (abs(low[indl]) + 1)
      hitl <- which(x[indl]-low[indl] <= toll)
      if (length(hitl) > 0) {
         flast <- f
         ipivot[indl[hitl]] <- -1
      }
   }
##---------------------------------------------------------------------
   indu <- which((ipivot == 0) & (p > 0));
   if (length(indu) > 0) {
      tolu <- 10 * eps * (abs( up[indu]) + 1)
      hitu <- which(up[indu]-x[indu]  <= tolu)
      if (length(hitu) > 0) {
         flast <- f
         ipivot[indu[hitu]] <- 1
      }
   }
##---------------------------------------------------------------------
   flast1  <- flast
   ipivot1 <- ipivot
   list(flast1 = flast, ipivot1=ipivot)
}

msolve <- function(g, upd1, ireset, first, d) {
##---------------------------------------------------------
## This routine acts as a preconditioning step for the
## linear conjugate-gradient routine.  It is also the
## method of computing the search direction from the
## gradient for the non-linear conjugate-gradient code.
## It represents a two-step self-scaled bfgs formula.
##---------------------------------------------------------
#%    cat("msolve, envjn$john=",envjn$john,"\n")
   envjn$john<-"msolve"
   if (upd1) {
#%       cat("upd1 is TRUE\n")
      y <- g/d
   } else {
      gsk <- as.numeric(crossprod(g, envjn$sk))
      if (ireset) {
         envjn$hg <- g/d 
         if (first) {
	    envjn$hyk   <- envjn$yk/d 
            ykhyk <- as.numeric(crossprod(envjn$yk, envjn$hyk))
         }
         ghyk <- as.numeric(crossprod(g, envjn$hyk)) 
         y    <- ssbfgs(envjn$sk,envjn$hg,envjn$hyk,envjn$yksk,ykhyk,gsk,ghyk)
      } else {
         envjn$hg <- g/d 
         if (first) {
            envjn$hyk   <- envjn$yk/d 
            envjn$hyr   <- envjn$yr/d 
            envjn$yksr  <- as.numeric(crossprod(envjn$yk, envjn$sr))
            ykhyr <- as.numeric(crossprod(envjn$yk, envjn$hyr))
         }
         gsr <- as.numeric(crossprod(g, envjn$sr))
         ghyr <- as.numeric(crossprod(g, envjn$hyr))
         if (first) {
            yrhyr <- as.numeric(crossprod(envjn$yr, envjn$hyr))
         }
         envjn$hg <- ssbfgs(envjn$sr,envjn$hg,envjn$hyr,envjn$yrsr,yrhyr,gsr,ghyr)
         if (first) {
            envjn$hyk <- ssbfgs(envjn$sr,envjn$hyk,envjn$hyr,envjn$yrsr,yrhyr,envjn$yksr,ykhyr)
         }
         ykhyk <- as.numeric(crossprod(envjn$hyk, envjn$yk))
         ghyk  <- as.numeric(crossprod(envjn$hyk, g)) 
         y     <- ssbfgs(envjn$sk,envjn$hg,envjn$hyk,envjn$yksk,ykhyk,gsk,ghyk) 
      }
   }
   amsolve<-list(y=y) ## returns y
}

ndia3<-function(e, v, gv, r, vgv){
##---------------------------------------------------------
## update the preconditioning matrix based on a diagonal 
## version of the bfgs quasi-newton update.
##---------------------------------------------------------
  tol    <- 1e-6
  vr     <- as.numeric(crossprod(v,r))
#%   cat("ndia3: vgv, vr, then e:",vgv, vr,"\n")
#%   print(e)
  if (abs(vr)>tol && abs(vgv)>tol) {
    e <- e - (r*r)/vr + (gv*gv)/vgv # CAUTION! May be crossprod
    ind  <- which(e < tol)
    if (length(ind)>0){
       e[ind] <- 1
    }
  }
#%   print(e)
  e # Need to return the object? Possibly not!!?
}

norm2 <- function (x) {
##---------------------------------------------------------
## compute the 2 norm of vector x
##---------------------------------------------------------
   res<-sqrt(as.numeric(crossprod(x)))
}
ssbfgs <- function (s, hv, hy, ys, yhy, vs, vhy){
##---------------------------------------------------------
## self-scaled bfgs quasi-Newton update 
## (used by preconditioner)
##---------------------------------------------------------
  delta <- (1 + yhy/ys)*vs/ys - vhy/ys
  beta  <- -vs/ys
  z     <- hv + delta*s + beta*hy
}

step1 <- function(f, gtp, smax) {
##---------------------------------------------------------
## step1 returns the length of the initial step to be 
## taken along the vector p in the next linear search.
##---------------------------------------------------------
## [fm is supposed to be an estimate of the optimal function value]
   eps<-.Machine$double.eps
   fm <- 0
   d  <- abs(f-fm)
   alpha <- min(1, smax)
   if ((2*d <= (-gtp)) && (d >= eps)) {
      alpha = min(-2*d/gtp, smax)
   }
   alpha # need to ensure returned value
}

stpmax <- function(stepmx, pe, x, p, ipivot, low, up) {
##------------------------------------------------
## compute the maximum allowable step length
## (spe is the standard (unconstrained) max step)
##------------------------------------------------
   ## cat("stpmax: stepmx, pe:",stepmx, pe,"\n")

   spe  <- stepmx / pe 
   al   <- spe 
   au   <- spe 
##------------------------------------------------
   indl <- which(ipivot==0 & p < 0) 
   if ( length(indl)>0) {
      tl   <- low[indl] - x[indl] 
      al   <- min(tl/p[indl]) 
   }
##------------------------------------------------
   indu <- which(ipivot==0 & p > 0) 
   if (length(indu)>0) {
      tu   <- up[indu] - x[indu] 
      au   <- min(tu/p[indu]) 
   }
##------------------------------------------------
   spe  <- min(c(spe, al, au)) 
}

ztime<-function(x,ipivot) {
   ## ---------------------------------------------------------
   ##  this routine multiplies the vector x by the 
   ##  constraint matrix z
   ## ---------------------------------------------------------
   ind <- which(ipivot!=0);
   x[ind] <- 0
   x1 <- x;
}
