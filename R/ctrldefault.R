##################################################################
ctrldefault <- function(npar) { 
# THIS IS FULL VERSION FOR optimx
#
# allmeth includes methods that are experimental or superceded
      allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3c", "Rcgmin", "Rtnmin", "Rvmmin", "snewton", "snewtonm",
                 "spg", "ucminf", "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb", 
                 "hjn", "lbfgs", "subplex", "ncg", "nvm", "mla", 
                 "slsqp", "tnewt", "anms", "pracmanm", "nlnm", "snewtm")

      truename <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3c", "Rcgmin", "Rtnmin", "Rvmmin", "snewton", "snewtonm",
                 "spg", "ucminf", "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb", 
                 "hjn", "lbfgs", "subplex", "ncg", "nvm", "mla", 
                 "slsqp", "tnewton", "anms", "nelder_mead", "neldermead", "snewtm")

#  allpkg has package where element of allmeth is found
      allpkg <-  c("stats", "stats", "stats", "stats", "stats", "stats",
                "lbfgsb3c", "optimx", "optimx", "optimx", "optimx", "optimx",
                "BB", "ucminf", "minqa", "minqa", "minqa", "dfoptim", "dfoptim", 
                 "optimx", "lbfgs", "subplex", "optimx", "optimx", "marqLevAlg", 
                 "nloptr", "nloptr", "pracma", "pracma", "nloptr", "optimx")

## These are DEFAULTS. They may be nonsense in some contexts.

      mostmeth <- c("Nelder-Mead", "nlm", "nlminb", 
                "lbfgsb3c", "Rtnmin", "snewtm",
                 "spg", "ucminf", "bobyqa", "nmkb", 
                 "subplex", "ncg", "Rcgmin", "nvm", 
                 "Rvmmin", "mla", "slsqp", "tnewt",
                 "pracmanm", "nlnm")

      nogrmeth <- c("Nelder-Mead", "newuoa", "bobyqa", "uobyqa", "nmkb", "hjkb",
                 "hjn", "subplex", "anms", "pracmanm", "nlnm")


      grmeth <- c("BFGS", "CG", "L-BFGS-B", "nlm", "nlminb", 
                "lbfgsb3c", "Rcgmin", "Rtnmin", "Rvmmin", 
                 "spg", "ucminf", "lbfgs", "ncg", "nvm", "mla", 
                 "slsqp", "tnewt")

     hessmeth <- c("nlm", "nlminb", "snewton", "snewtonm", "snewtm")
# NOTE: snewtm and snewtonm are synonyms



#  allpkg has package where element of allmeth is found
#      mostpkg <-  c("stats", "stats", "stats",
#                "lbfgsb3c", "optimx", "optimx", 
#                "BB", "ucminf", "minqa", "dfoptim",  
#                "subplex", "optimx", "optimx", "marqLevAlg", 
#                "nloptr", "pracma", "nloptr")

# Adjust for packages not installed. 
# !! COMMENTED OUT TO AVOID UNNECESSARY WORK
     #  uap <- unique(allpkg)
     #  nu <- length(uap)
     #  tu <- rep(NA, nu) # to define array for indicator of installed packages
     #  for (i in 1:nu){ tu[i] <- require(uap[i], character.only=TRUE, quietly=TRUE) }
     #  if (all(tu) ) {
     #    # Change allmeth to xxxmeth, allpkg to xxxpkg if experimental pkgs wanted
          OKmeth <- allmeth
          OKpkg <- allpkg
     #  } 
     #  else {
     #    badp <- uap[ - which(tu) ]
     #    for (i in 1:length(badp)) { warning("Package ",badp[i]," not installed") }
     #    uap <- uap[ which(tu) ]
     #    OK <- (allpkg %in% uap)
     #    badm <- which( ! OK )
     #    OKmeth <- allmeth[ - badm ] # leave only packages not 
     #    OKpkg <- allpkg[ OK]
     # }
      weakmeth <- c("snewton", "uobyqa")

      bdmeth <- c("L-BFGS-B", "nlminb", "lbfgsb3c", "Rcgmin", "Rtnmin", "nvm",  
                "Rvmmin", "bobyqa", "nmkb", "hjkb", "hjn", "snewtonm", "ncg", 
                "slsqp", "tnewt", "nlnm", "snewtm", "spg")
                 # snewtonmb added 20220210, removed 20230625 for snewtm

      bdmeth <- bdmeth[ which(bdmeth %in% OKmeth) ]
   

      maskmeth <- c("Rcgmin", "nvm", "hjn", "ncg", "snewtonm", "nlminb", "L-BFGS-B") 
      maskmeth <- maskmeth[ which(maskmeth %in% OKmeth) ]
  
#     valid gradient approximations
      grapprox <- c("grfwd", "grcentral", "grback", "grnd", "grpracma") 
 
# offset changed from 100 to 1000 on 180710
      ctrl.default <- list(
        acctol = 0.0001, # used for acceptable point test in backtrack linesearch
        all.methods = FALSE, # we do NOT want all methods to be the default
        allmeth = OKmeth, # to define the set of all methods
        truename = truename, # true names
        allpkg = OKpkg, # to list all the packages required
        appgr = FALSE, # assume we are NOT using a numerical approximation to the gradient
#        avoidmeth = avoidmeth, # methods to avoid using
      	badval = (0.5)*.Machine$double.xmax, # use this value as a flag of non-computable item
        bdmeth = bdmeth, # list of bounds-constrained methods
        bigval = .Machine$double.xmax*0.01, # a very large number (note smaller than badval)
        checkgrad = FALSE, # Rcgmin and Rvmmin check analytic gradient
        defgrapprox = "grfwd", # use forward approximation as default. Could argue for grcentral
        defmethod = "Nelder-Mead", # same default as optim(). Do we want this?
        defstep=1, # default stepsize is 1 (Newton stepsize)
        dowarn = TRUE,  # generally want to turn on warnings
        eps = 1e-07, # a tolerance for a small quantity (about single precision level)
        epstol = .Machine$double.eps, # but this is the machine epsilon
        fnscale = 1.0, # Normally scale function by 1. -1 maximizes functions
        follow.on = FALSE, # for optimx when multiple methods
        grapprox = grapprox, # set of valid gradient approximation functions
        grtesttol=(.Machine$double.eps)^(1/3), # a test tolerance for equality for numeric and 
        # analytic gradients
        have.bounds = FALSE, # normally have UNCONSTRAINED function
        have.masks = FALSE, # normally no masks
        hessasymtol = 0.0001, # tolerance for testing Hessian asymmetry. Note that it 
        # cannot be too small or we get too many false positives
        hesspkg="numDeriv", # default package for origin of hessian()
        hesstesttol=(.Machine$double.eps)^(1/3), # See grtesttol. This is for Hessian approx.
        jacpkg="numDeriv", # default package for origin of jacobian()
        keepinputpar = FALSE, # When TRUE do NOT allow bounds check to change parameter values
        kkt = TRUE, # Normally test KKT conditions. May take a LONG time.
        kkttol = 0.001, # tolerance for testing KKT small gradient condition 
        kkt2tol = 1.0E-6, # tolerance for testing KKT curvature condition
        maskmeth = maskmeth, # list of methods that allow masks (fixed parameters)
        maximize = FALSE, # normally MINIMIZE (see fnscale)
        maxit = 500*round(sqrt(npar+1)), # limit on number of iterations or gradient evaluations
        maxfeval = 5000*round(sqrt(npar+1)), # limit on function evaluations
        mostmeth = mostmeth,
        parchanged = FALSE, # set TRUE when bounds check has changed parameter values
        # ?? do we want this as a CONTROL? It is returned how??
        parscale = rep(1, npar), # vector of scaling factors for parameters. Try to get
        # scaled parameters to have magnitude in range (1, 10)
        reltest = 100.0, # used for equality test (a + offset) == (b + offset)
        save.failures = TRUE, # ?? where used. optimx() saves failed runs. opm?
      	scaletol = 3, 
        starttests = FALSE, 
        stepinc = 10.0, # used for lambda increase etc.
        steplen0 = 0.75, 
        stepmax = 5,
        stepmin = 0,
        stepredn = 0.2, # also used for lambda decrease
        stopbadupdate = FALSE,
        stoponerror = FALSE, # Don't stop with try-error normally
        tol = 0, 
        trace = 0,
        watch = FALSE,
        weakmeth = weakmeth # methods that should not be used
      )
      ctrl.default
}
##################################################################

dispdefault <- function(ctrl) { # Display the control vector using cat and print
  cat("Control vector (ctrl) multiple element items first:\n")
  cat("mostmeth:")
  print(ctrl$mostmeth)
  cat("allmeth:")
  print(ctrl$allmeth)
  cat("allpkg:")
  print(ctrl$allpkg)
  cat("bdmeth:")
  print(ctrl$bdmeth)
  cat("maskmeth:")
  print(ctrl$maskmeth)
  cat("mostmeth:")
  print(ctrl$mostmeth)
  cat("parscale:")
  print(ctrl$parscale)
  cat("acctol=",ctrl$acctol,"  all.methods=",ctrl$all.methods,"  badval=",ctrl$badval,
    "  bigval=",ctrl$bigval,"\n")
  cat("defgrapprox=",ctrl$defgrapprox,"  defmethod=",ctrl$defmethod,"  defstep=",ctrl$defstep,"\n")
  cat("dowarn=",ctrl$dowarn,"  eps=",ctrl$eps,"  epstol=",ctrl$epstol,"  fnscale=",ctrl$fnscale,"\n")
  cat("grtesttol=",ctrl$grtesttol,"  have.bounds=",ctrl$have.bounds,"  hessasymptol=",ctrl$hessasymptol,"\n")
  cat("hesttesttol=",ctrl$hesstesttol,"  keepinputpar=",ctrl$keepinputpar,"  kkt=",ctrl$kkt,"\n")
  cat("kkttol=",ctrl$kkttol,"  kkt2tol=",ctrl$kkt2tol,"  maximize=",ctrl$maximize,"\n")
  cat("maxit=",ctrl$maxit,"  maxfeval=",ctrl$maxfeval,"  offset=",ctrl$offset,
     "  parchanged=",ctrl$parchanged,"\n")
  cat("reltest=",ctrl$reltest,"  save.failures=",ctrl$save.failures,"  scaletol=",ctrl$scaletol,"\n")
  cat("  stepinc=",ctrl$stepinc,"  steplen0=",ctrl$steplen0,"  stepmax=",ctrl$stepmax,
     "  stepmin=",ctrl$stepmin,"stepredn=",ctrl$stepredn,"\n")
  cat("  stopbadupdate=",ctrl$stopbadupdate,"  tol=",ctrl$tol,"\n")
  cat("trace=",ctrl$trace,"  watch=",ctrl$watch,"\n")  
  cat("=========================================\n")  
  return(0)
}

