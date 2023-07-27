bmchk <- function(par, lower = NULL, upper = NULL, 
    bdmsk = NULL, trace = 0, offset = 100.0, shift2bound = TRUE) {
## 
## 2022-3-14 Issues:
## - bdmsk on input should be checked SEPARATELY from bounds and masks for consistency
## - for setting mask (maskadded flag) want to use offset (formerly tol method)
## 
    ## Bounds and masks check
    #  !!? check use of par and bvec -- can we simplify?
    # 20101031 -- issue of bounds not working correctly
    #  - -Inf seems to upset bounds
    #  - no upper bounds gives troubles (applies to Rcgmin too!)
    #
    # Input:
    #  par = a vector containing the starting point
    #  lower = vector of lower bounds on parameters
    #  upper = vector of upper bounds on parameters
    # Note: free parameters outside bounds will be adjusted to
    #   bounds.
    # bdmsk = control vector for bounds and masks. Parameters
    #   for which bdmsk are 1
    # are unconstrained or 'free', those with bdmsk 0 are
    #   masked i.e., fixed.
    # For historical reasons, we use the same array as an
    #   indicator that a parameter is at a lower bound (-3)
    #   or upper bound (-1)
    # trace = control of output: 0 for none (default), >0 for
    #   output
    # offset = shift value to test for closeness of upper and lower bounds
    #   (a + offset) =?= (b + offset) avoids use of a tolerance
    # shift2bound = TRUE if we adjust par values so they are
    #   feasible
    ##
    # Output:
    #    A list with components:
    #     bvec: The parameters adjusted to the nearest bound.
    #     bdmsk: adjusted input masks
    #     bchar: indicator for humans -- "-","L","F","U","+","M","?","!"
    #        for out-of-bounds-low, lower bound, free, 
    #            upper bound, out-of-bounds-high, masked (fixed),
    #            unknown (?), inadmissible (!)
    #     lower: adjusted lower bounds
    #     upper: adjusted upper bounds
    #     nolower: TRUE if no lower bounds, FALSE otherwise
    #     noupper: TRUE if no upper bounds, FALSE otherwise
    #     bounds:  TRUE if any bounds, FALSE otherwise ?? at moment confused about masks being bounds
    #     admissible: TRUE if admissible, FALSE if not
    #        No lower bound exceeds an upper bound. That is the 
    #        bounds themselves are sensible. This condition has 
    #        nothing to do with the starting parameters.
    #     maskadded: TRUE if a mask is added, FALSE if not
    #        This implies very close upper and lower bounds for 
    #        the parameters. See the code for the implementation.
    #     parchanged: TRUE if parameters changed, FALSE if not
    #        parchanged = TRUE means that parameters are 
    #        INFEASIBLE, or they would not be changed.
    #     onbound: TRUE if any parameter equal to a bound
    #     feasible: TRUE if feasible, FALSE otherwise
    #
    ########## length of vectors #########
    n <- length(par)
    bvec <- par
    ############# bounds and masks ################
    # set default masks if not defined
    bchar <- rep(" ",n) # initialize indicators to "blank"
    if (is.null(bdmsk)) { bdmsk <- rep(1, n) } # initialize to free if null
    if (trace > 2) { cat("bdmsk:"); print(bdmsk) }
    # check if there are bounds (non null and no finite values)
    if (is.null(lower) || !any(is.finite(lower))) nolower <- TRUE
    else nolower <- FALSE
    if (is.null(upper) || !any(is.finite(upper))) noupper <- TRUE
    else noupper <- FALSE
    if (nolower && noupper && all(bdmsk == 1)) bounds <- FALSE
    else bounds <- TRUE  # bounds indicator now set
    if (trace > 2) cat("Bounds: nolower = ", nolower, "  noupper = ", noupper, 
            " bounds = ", bounds, "\n")
    if (nolower) lower <- rep(-Inf, n)
    if (noupper) upper <- rep(Inf, n)

## adjust tolerance for masks and parameters ON bounds
#    if (is.null(tol) || tol <= 0.0) {
#       tol <- .Machine$double.eps * max(abs(par), 1) # use par as bounds Inf
#    }
    if (trace > 1) {
        cat("Initial parameters:"); print(par)
    }
    ######## check bounds and masks #############
    parchanged <- FALSE  # must be set BEFORE if (bounds) ...; 
                         # initialized to indicate parameters NOT changed
    feasible <- TRUE    # initially assume parameters are feasible
    admissible <- TRUE  # similarly admissible (must set before we look at bounds)
    maskadded <- FALSE  # indicate no masks added
    onbound <- FALSE    # indicate NOT on bounds
    if (bounds) {  # Make sure to expand lower and upper
        if (!nolower & (length(lower) < n)) 
        {   if (length(lower) == 1) { lower <- rep(lower, n) }
            else { stop("1<length(lower)<n") }
        }  # else lower OK
        if (!noupper & (length(upper) < n)) 
        {   if (length(upper) == 1) { upper <- rep(upper, n) }
            else { stop("1<length(upper)<n") }
        }  # else upper OK
        # At this point, we have full bounds in play
        ######## check admissibility ########
        if (any(lower[which(bdmsk != 0)] > upper[which(bdmsk != 0)])) admissible <- FALSE
        if (trace > 0) cat("admissible = ", admissible, "\n")
        if ( any((upper+offset) == (lower + offset)) ) { # essentially masked
            makemask<-which((upper + offset) == (lower + offset))
            if (trace > 0) {
               cat("Imposing mask as lower ~= upper for following parameters\n")
               print(makemask)
            }
            # force parameters to the masked values
            bvec[makemask] <- lower[makemask] + 0.5*(upper[makemask] - lower[makemask])
            if (any(bvec != par)) { 
               parchanged <- TRUE 
               warning("Masks (fixed parameters) set by bmchk due to tight bounds. CAUTION!!")
            }
            bdmsk[makemask] <- 0 # set bchar below
            maskadded <- TRUE
        }
        if (trace > 0) cat("maskadded = ", maskadded, "\n")
        ######## check feasibility ########
        if (admissible) { # This implementation as a loop, but try later to vectorize
            for (i in 1:n) {
                if (bdmsk[i] == 0) {
                  bchar[i] <- "M" # NOTE: do not change masked parameters, 
                                  # even if out of bounds
                  if (!nolower) { # there is a lower bound
                    if (bvec[i] < lower[i]) {
                      if (trace > 0) { 
                        cat("WARNING: ", bvec[i], " = MASKED x[", 
                          i, "] < lower bound = ", lower[i], "\n")
                      }
                      feasible <- FALSE
                    }
                  }
                  if (!noupper) { # there is an upper bound
                    if (bvec[i] > upper[i]) {
                      if (trace > 0){
                        cat("WARNING: ", bvec[i], " = MASKED x[", 
                          i, "] > upper bound = ", upper[i], "\n")
                      }
                      feasible<-FALSE
                    }
                  }
                }
                else { # par[i] not masked, so must be free or active constraint
                  if (!nolower) { # there are lower bounds
                    if (bvec[i] <= lower[i]) { # on or below lower bound
                      if ((offset+bvec[i]) < (offset+lower[i])){
                         bdmsk[i] <- -3.5 # OUT OF BOUNDS LOW
                         bchar[i] <- "-"
                         feasible<-FALSE
                         if (shift2bound) { # allow shift of parameter
                            parchanged <- TRUE
                            if (trace > 0) cat("WARNING: x[", i, "], set ", 
                                 bvec[i], " to lower bound = ", lower[i], "\n")
                            bvec[i] <- lower[i]
                            bdmsk[i] <- -3 
                            bchar[i] <- "L"
                         }
                      } else {
                         bdmsk[i] <- -3  # active lower bound
                         bchar[i] <- "L"
                      }
                    }
                  } # ! nolower
                  if (!noupper) { # there are upper bounds
                    if (bvec[i] >= upper[i]){ # on or above upper bound
                      if ((offset+bvec[i]) > (offset+upper[i])) {
                         bdmsk[i] <- -0.5 # OUT OF BOUNDS HIGH
                         bchar[i] <- "+"
                         feasible<-FALSE
                         if (shift2bound) {
                           parchanged <- TRUE
                           if (trace > 0) cat("WARNING: x[", i, "], set ",
                                  bvec[i], " to upper bound = ", upper[i], "\n")
                           bvec[i] <- upper[i]
                           bdmsk[i] <- -1   # active upper bound
                           bchar[i] <- "U"
                         }
                      } else { 
                         bdmsk[i] <- -1  # active upper bound
                         bchar[i] <- "U"
                      }
                    }
                  } # noupper
               }  # end not masked
               ## Should NOT proceed with optimization if feasible and parchanged 
               ## both FALSE (out of bounds)
            }  # end loop for bound/mask check
        } # if admissible
        else { # NOT admissible. 
           bchar <- rep("!",n)
        }
        if (trace > 1) {
           cat("lower:"); print(lower)
           cat("upper:"); print(upper)
        }
    } # if bounds
    if (trace > 0) 
        cat("parchanged = ", parchanged, "\n")
    if ((trace > 1) && parchanged) {
       cat("changed par:"); print(par)
    }
    if (any((bvec+offset) == (lower+offset)) || any((bvec+offset) == (upper+offset)) ) {
        onbound <- TRUE
        if (trace > 0) cat("At least one parameter is on a bound\n")
    }
    if ((trace > 0) && (! is.null(attr(bvec,"status"))) ) {
        cat("existing parameter status:"); print(attr(bvec,"status")) 
    }
    attr(bvec,"status")<-bchar # (re)set status
    ############## end bounds check #############
    bcout <- list(bvec, bdmsk, bchar, lower, upper, nolower, noupper, 
        bounds, admissible, maskadded, parchanged, feasible, onbound)
    names(bcout) <- c("bvec", "bdmsk", "bchar", "lower", "upper", "nolower", "noupper", 
       "bounds", "admissible", "maskadded", "parchanged", "feasible", "onbound")
    # Note bdmsk, lower and upper are returned because they are modified (length, etc.)
    return(bcout)
}  ## end of bmchk.R
