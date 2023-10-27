ncg <- function(par, fn, gr, bds, control = list()) {
    # control defaults -- idea from spg
    ctrl <- list(maxit = 500, trace = 0, eps = 1e-07, dowarn = TRUE, tol=0, stepredn=0.15)
    namc <- names(control)
#    if (!all(namc %in% names(ctrl))) 
#        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control[namc]
    npar<-length(par)
    if (ctrl$tol == 0) tol <- npar * (npar * .Machine$double.eps)  # for gradient test.  
      # Note -- possible integer overflow if n*n*d.eps
    else tol<-ctrl$tol
    maxit <- ctrl$maxit  # limit on function evaluations
    trace <- ctrl$trace  # 0 for no output, >0 for output (bigger => more output)
    if (trace > 2) cat("trace = ", trace, "\n")
    eps <- ctrl$eps
    grNULL <- is.null(gr)
    dowarn <- ctrl$dowarn  #
    # gr MUST be provided
    if (is.null(gr)) { 
      stop("A gradient calculation (analytic or numerical) MUST be provided for ncg") 
    } else { mygr<-gr }
    ## Set working parameters (See CNM Alg 22)
    if (trace > 0) {
        cat("ncg -- J C Nash 2023 - bounds constraint version of new CG\n")
        cat("an R implementation of Alg 22 with Yuan/Dai modification\n")
    }
    bvec <- par  # copy the parameter vector
    maxfeval <- round(sqrt(npar + 1) * maxit)  # change 091219
    ig <- 0  # count gradient evaluations
    ifn <- 1  # count function evaluations (we always make 1 try below)
    stepredn <- ctrl$stepredn
    if (trace > 0) cat("stepredn =",stepredn,"\n") 
    acctol <- 1e-04  # acceptable point tolerance
    reltest <- 100  # relative equality test
    ceps <- .Machine$double.eps * reltest
    accpoint <- as.logical(FALSE)  # so far do not have an acceptable point
    cyclimit <- min(2.5 * npar, 10 + sqrt(npar))  #!! upper bound on when we restart CG cycle
    # set default masks if not defined
    if (is.null(bds)) {
        bdmsk <- rep(1, npar)
    }
    else bdmsk <- bds$bdmsk
    bounds <- bds$bounds
    if (trace > 2) { cat("bdmsk:"); print(bdmsk)}
    if (trace > 2) cat("Bounds: nolower = ", bds$nolower, "  noupper = ", bds$noupper, 
            " bounds = ", bounds, "\n")
    # Initial function value -- may NOT be at initial point specified by user.
    if (trace > 2) {cat("Try function at initial point:");  print(bvec)  }
    f <- try(fn(bvec), silent = FALSE)  # Compute the function at initial point.
    if (trace > 0) { cat("Initial function value=", f, "\n") }
    if (inherits(f,"try-error")) {
        msg <- "Initial point is infeasible."
        if (trace > 0) 
            cat(msg, "\n")
        ans <- list(par, NA, c(ifn, 0), 2, msg, bdmsk)
        names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
        return(ans)
    }
    fmin <- f
    if (trace > 0) cat("Initial fn=", f, "\n")
    if (trace > 2) print(bvec)
    # Start the minimization process
    keepgoing <- TRUE
    msg <- "not finished"  # in case we exit somehow
    oldstep <- 0.8  #!! 2/3 #!!?? WHY?
    fdiff <- NA  # initially no decrease
    cycle <- 0  # !! cycle loop counter
    ####################################################################
    while (keepgoing) { # main loop -- must remember to break out of it!!
        t <- as.vector(rep(0, npar))  # zero step vector
        c <- t  # zero 'last' gradient
        while (keepgoing && (cycle < cyclimit)) {    ## cycle loop
            cycle <- cycle + 1
            if (trace > 0) 
                cat(ifn, " ", ig, " ", cycle, " ", fmin, "  last decrease=", 
                  fdiff, "\n")
            if (trace > 2) { print(bvec); cat("\n") }
            if (ifn > maxfeval) {
                msg <- paste("Too many function evaluations (> ", 
                  maxfeval, ") ", sep = "")
                if (trace > 0) 
                  cat(msg, "\n")
                ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  # 1 indicates not converged in function limit
                names(ans) <- c("par", "value", "counts", "convergence", 
                  "message", "bdmsk")
                return(ans)
            }
            par <- bvec  # save best parameters
            ig <- ig + 1
            if (ig > maxit) {
                msg <- paste("Too many gradient evaluations (> ", maxit, ") ", sep = "")
                if (trace > 0) cat(msg, "\n")
                ans <- list(par, fmin, c(ifn, ig), 1, msg, bdmsk)  
                       # 1 indicates not converged in function or gradient limit
                names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
                return(ans)
            }
            g <- mygr(bvec)
            if (bounds) { ## Bounds and masks adjustment of gradient ##
                  ## first try with looping -- later try to vectorize
                  if (trace > 2) { cat("bdmsk:"); print(bdmsk) }
                  for (i in 1:npar) {
                    if ((bdmsk[i] == 0)) {  g[i] <- 0 }# masked: gradient component is zero
                    else {
                      if (bdmsk[i] == 1) {
                        if (trace > 1) cat("Parameter ", i, " is free\n")
                      }
                      else {
                        if ((bdmsk[i] + 2) * g[i] < 0) {
                          # test for -ve gradient at upper bound, +ve at lower bound
                          g[i] <- 0  # in which case active mask or constraint and zero gradient component
                        }
                        else {
                          bdmsk[i] <- 1  # freeing parameter i
                          if (trace > 1) cat("freeing parameter ", i, "\n")
                        }
                      }
                    }
                  }  # end masking loop on i
                  if (trace > 2) { 
                    cat("bdmsk adj:\n"); print(bdmsk)
                    cat("proj-g:\n") ; print(g)
                  }
                }  # end if bounds
            ## end bounds and masks adjustment of gradient
            g1 <- sum(g * (g - c))  # gradient * grad-difference
            g2 <- sum(t * (g - c))  # oldsearch * grad-difference
            gradsqr <- sum(g * g)
            if (trace > 1) {
                cat("Gradsqr = ", gradsqr, " g1, g2 ", g1, " ", 
                  g2, " fmin=", fmin, "\n")
            }
            c <- g  # save last gradient
            g3 <- 1  # !! Default to 1 to ensure it is defined -- t==0 on first cycle
            if (gradsqr > tol * (abs(fmin) + reltest)) {
                if (g2 > 0) {
                  betaDY <- gradsqr/g2
                  betaHS <- g1/g2
                  g3 <- max(0, min(betaHS, betaDY))  # g3 is our new 'beta' !! Dai/Yuan 2001, (4.2)
                }
            }
            else {
                msg <- paste("Very small gradient -- gradsqr =", gradsqr, sep = " ")
                if (trace > 0) cat(msg, "\n")
                keepgoing <- FALSE  # done loops -- should we break ??
                break  # to leave inner loop
            }
            if (trace > 2) cat("Betak = g3 = ", g3, "\n")
            if (g3 == 0 || cycle >= cyclimit) {
                # we are resetting to gradient in this case
                if (trace > 0) {
                  if (cycle < cyclimit) cat("Yuan/Dai cycle reset\n")
                  else cat("Cycle limit reached -- reset\n")
                }
                fdiff <- NA
                cycle <- 0
                break  #!!
            }
            else {  # drop through if not Yuan/Dai cycle reset
                t <- t * g3 - g  # t starts at zero, later is step vector
                gradproj <- sum(t * g)  # gradient projection
                if (trace > 1) cat("Gradproj =", gradproj, "\n")
                if (bounds) {  ## Adjust search direction for masks
                    if (trace > 2) {cat("t:\n"); print(t)}
                    t[which(bdmsk <= 0)] <- 0  # apply mask constraint
                    if (trace > 2) { cat("adj-t:\n");  print(t)}
                    ## end adjust search direction for masks
                  }  # end if bounds
                # ?? Why do we not check gradproj size??
                ########################################################
                ####                  Line search                   ####
                OKpoint <- FALSE
                if (trace > 2) cat("Start linesearch with oldstep=", oldstep, "\n")
                steplength <- oldstep * 1.5  #!! try a bit bigger
                f <- fmin
                changed <- TRUE  # Need to set so loop will start
                while ((f >= fmin) && changed) {
                  if (bounds) { # Box constraint -- adjust step length
                      for (i in 1:npar) { # loop on parameters -- vectorize??
                        if ((bdmsk[i] == 1) && (t[i] != 0)) {
                            # only concerned with free parameters and non-zero search dimension
                            if (t[i] < 0) { # going down. Look at lower bound
                              trystep <- (bds$lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
                            }
                            else {  # going up, check upper bound
                              trystep <- (bds$upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
                            }
                            if (trace > 2) cat("steplength, trystep:", steplength, trystep, "\n")
                            steplength <- min(steplength, trystep)  # reduce as necessary
                          }  # end steplength reduction
                      }  # end loop on i to reduce step length
                      if (trace > 1) cat("reset steplegth=", steplength, "\n")
                      # end box constraint adjustment of step length
                    }  # end if bounds
                  bvec <- par + steplength * t
                  changed <- (!identical((bvec + reltest), (par + reltest)))
                  if (changed) { # compute newstep, if possible
                    f <- fn(bvec)  # Because we need the value for linesearch, don't use try()
                    # instead preferring to fail out, which will hopefully be unlikely.
                    ifn <- ifn + 1
                    if (is.na(f) || (!is.finite(f))) { warning("ncg - undefined function")
                      f <- .Machine$double.xmax
                    }
                    if (f < fmin) {
                      f1 <- f  # Hold onto value
                    }
                    else {
                      savestep <- steplength
                      steplength <- steplength * stepredn
#                      cat("stepredn:",stepredn,"\n")
#                      cat("savestep:",savestep,"\n")
#                      cat("steplength:",steplength,"\n")
                      if (steplength >=savestep) changed <- FALSE
                      if (trace > 0) cat("*")
                    }
                  }
                }  # end while
                changed1 <- changed  # Change in parameters occured in step reduction
                if (changed1) {  ## ?? should we check for reduction? or is this done in if (newstep >0) ?
                    newstep <- 2 * (f - fmin - gradproj * steplength)  # JN 081219 change
                    if (newstep > 0) {
                      newstep = -(gradproj * steplength * steplength/newstep)
                    } 
                    else newstep <- 2*steplength # 20220627 try a doubling
                    if (bounds) { # Box constraint -- adjust step length
                        for (i in 1:npar) {  # loop on parameters -- vectorize??
                          if ((bdmsk[i] == 1) && (t[i] != 0)) {
                              # only concerned with free parameters and non-zero search dimension
                              if (t[i] < 0) { # going down. Look at lower bound
                                trystep <- (bds$lower[i] - par[i])/t[i]  # t[i] < 0 so this is positive
                              }
                              else { # going up, check upper bound
                                trystep <- (bds$upper[i] - par[i])/t[i]  # t[i] > 0 so this is positive
                              }
                              if (trace > 2) cat("newstep, trystep:", newstep, trystep, "\n")
                              newstep <- min(newstep, trystep)  # reduce as necessary
                          }  # end newstep reduction
                        }  # end loop on i to reduce step length
                        if (trace > 2) cat("reset newstep=", newstep, "\n")
                        # end box constraint adjustment of step length
                    }  # end if bounds
                    bvec <- par + newstep * t
                    changed <- (!identical((bvec + reltest), (par + reltest)))
                    if (changed) {
                      f <- fn(bvec)
                      ifn <- ifn + 1
                    }
                    if (trace > 2) cat("fmin, f1, f: ", fmin, f1, f, "\n")
                    if (f < min(fmin, f1)) { # success
                      OKpoint <- TRUE
                      accpoint <- (f <= fmin + gradproj * newstep * 
                        acctol)
                      fdiff <- (fmin - f)  # check decrease
                      fmin <- f
                      oldstep <- newstep  # !! save it
                    }
                    else {
                      if (f1 < fmin) {
                        bvec <- par + steplength * t  # reset best point
                        accpoint <- (f1 <= fmin + gradproj * steplength * acctol)
                        OKpoint <- TRUE  # Because f1 < fmin
                        fdiff <- (fmin - f1)  # check decrease
                        fmin <- f1
                        oldstep <- steplength  #!! save it
                      }
                      else { # no reduction
                        fdiff <- NA
                        accpoint <- FALSE
                      }  # f1<?fmin
                    }  # f < min(f1, fmin)
                    if (trace > 1) cat("accpoint = ", accpoint, " OKpoint = ", OKpoint, "\n")
                    if (!accpoint) {
                      msg <- "No acceptable point -- exit loop"
                      if (trace > 0) cat("\n", msg, "\n")
                      keepgoing <- FALSE
                      break  #!!
                    }
                }  # changed1
                else { # not changed on step redn
                  if (cycle == 1) { msg <- " Converged -- no progress on new CG cycle"
                    if (trace > 0) cat("\n", msg, "\n")
                    keepgoing <- FALSE
                    break  #!!
                  }
                }  # end else
#               if (! keepgoing) break # 2023-10-27 to exit
            }  # end of test on Yuan/Dai condition
            #### End line search ####
            if (bounds) { ## Reactivate constraints?? -- should check for infinite bounds
                 for (i in 1:npar) {
                    if (bdmsk[i] == 1) { # only interested in free parameters
                      if (is.finite(bds$lower[i])) { # JN091020 -- need to use abs in case bounds negative
                        if ((bvec[i] - bds$lower[i]) < ceps * (abs(bds$lower[i]) + 1)) {
                          # are we near or lower than lower bd
                          if (trace > 2) cat("(re)activate lower bd ", i, " at ", bds$lower[i], "\n")
                          bdmsk[i] <- -3
                        }  # end lower bd reactivate
                      }
                      if (is.finite(bds$upper[i])) { # JN091020 -- need to use abs in case bounds negative??
                        if ((bds$upper[i] - bvec[i]) < ceps * (abs(bds$upper[i]) + 1)) {
                           # are we near or above upper bd
                          if (trace > 2) cat("(re)activate upper bd ", i, " at ", bds$upper[i], "\n")
                          bdmsk[i] <- -1
                        }  # end lower bd reactivate
                      }
                    }  # end test on free params
                 }  # end reactivate constraints
            }  # end if bounds
#            if (! keepgoing) break # 2023-10-27 to exit
        }  # end of inner loop (cycle)
        if (oldstep < acctol) {
            oldstep <- acctol
        }  #   steplength
        if (oldstep > 1) { oldstep <- 1 } # Force no bigger than 1
        if (trace > 1) cat("End inner loop, cycle =", cycle, "\n")
#        if (! keepgoing) break # 2023-10-27 to exit
    }  # end of outer loop
    msg <- "ncg seems to have converged"
    if (trace > 0) cat(msg, "\n")
    #  par: The best set of parameters found.
    #  value: The value of 'fn' corresponding to 'par'.
    #  counts: number of calls to 'fn' and 'gr' (2 elements)
    # convergence: An integer code. '0' indicates successful
    #   convergence.
    #  message: A character string or 'NULL'.
    ans <- list(par, fmin, c(ifn, ig), 0, msg, bdmsk)
    names(ans) <- c("par", "value", "counts", "convergence", "message", "bdmsk")
    return(ans)
}  ## end of ncg version of ncg 230620

