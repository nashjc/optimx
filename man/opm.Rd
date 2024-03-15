\name{opm}
\alias{opm}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls multiple other
  R tools for optimization, including the existing optim() function tools.

  Because SANN does not return a meaningful convergence code
  (conv), \code{opm()} does not call the SANN method, but it can be invoked 
  in \code{optimr()}.

  There is a pseudo-method "ALL" that runs all available methods. Note that
  this is upper-case. This function is a replacement for optimx() from the
  optimx package.  \code{opm()} calls the \code{optimr()} function for each 
  solver in the \code{method} list.
}
\usage{
opm(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with a first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, whatever default actions are supplied by the methods
   specified will be used. However, some methods REQUIRE a gradient function, so 
   will fail in this case. \code{opm()} will generally return with \code{convergence}
   set to 9998 for such methods.

   If 'gr' is a character string, this character string will be taken to be the name
   of an available gradient approximation function. Examples are "grfwd", "grback", 
   "grcentral" and "grnd", with the last name referring to the default method of 
   package \code{numDeriv}.}
 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints. These are vectors.}
 \item{method}{A vector of the methods to be used, each as a character string.
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}. If \code{method}
       has just one element, \code{"ALL"} (capitalized), all available and 
       appropriate methods will be tried.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim().}
 \item{control}{A list of control parameters. See \sQuote{Details}. There is
    a spreadsheet \code{/inst/doc/optcontrol.xls} that is an ongoing attempt to
    document the different controls used in the various methods.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  For details of how \code{opm()} calls the methods, see the documentation
  and code for \code{optimr()}. The documentation and code for individual
  methods may also be useful. Note that some simplification of the calls
  may have been necessary, for example, to provide reasonable default values
  for method controls that are consistent across several methods, though this
  is not always possible to guarantee. The documentation for \code{optimr} and the source
  code of the quite simple routine \code{ctrldefault.R} may be useful.

  Some of the commonly useful elements of the \code{control} list are: 
  
  \describe{
    \item{\code{trace}}{Non-negative integer. If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing. trace = 0 gives no output 
      (To understand exactly what these do see the source code: higher 
      levels give more detail.)}
    \item{\code{maxfeval}}{For methods that can use this control, a limit
      on the number of function evaluations. This control is simply passed
      through. It is not checked by \code{opm}.}
    \item{\code{maxit}}{For methods that can use this control, a limit
      on the number of gradient evaluations or major iterations.}
    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} and \code{gr} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}. For methods from the set in
      \code{optim()}. Note potential conflicts with the control \code{maximize}.}
    \item{\code{parscale}}{A vector of scaling values for the parameters.
	Optimization is performed on \code{par/parscale} and these should be
	comparable in the sense that a unit change in any element produces
	about a unit change in the scaled value.For \code{optim}.}
    \item{\code{save.failures}}{ = TRUE (default) if we wish to keep "answers" from runs 
      where the method does not return convergence==0. FALSE otherwise.}
    \item{\code{maximize}}{ = TRUE if we want to maximize rather than minimize 
      a function. (Default FALSE). Methods nlm, nlminb, ucminf cannot maximize a
      function, so the user must explicitly minimize and carry out the adjustment
      externally. However, there is a check to avoid
      usage of these codes when maximize is TRUE. See \code{fnscale} below for 
      the method used in \code{optim} that we deprecate.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.  This is equivalent to setting \code{method="ALL"}}
    \item{\code{kkt}}{=FALSE if we do NOT want to test the Kuhn, Karush, Tucker
      optimality conditions. The default is generally TRUE. However, because the Hessian
      computation may be very slow, we set \code{kkt} to be FALSE if there are 
      more than than 50 parameters when the gradient function \code{gr} is not 
      provided, and more than 500
      parameters when such a function is specified. We return logical values \code{KKT1}
      and \code{KKT2} TRUE if first and second order conditions are satisfied approximately.
      Note, however, that the tests are sensitive to scaling, and users may need
      to perform additional verification. If \code{hessian} is TRUE, this overrides 
      control \code{kkt}.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkttol}}{= value to use to check for small gradient and negative
      Hessian eigenvalues. Default = .Machine$double.eps^(1/3) }
    \item{\code{kkt2tol}}{= Tolerance for eigenvalue ratio in KKT test of positive 
      definite Hessian. Default same as for kkttol }
    \item{\code{dowarn}}{= FALSE if we want to suppress warnings generated by \code{opm()} or
      \code{optimr()}. Default is TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.
      The value is then a signal of failure when execution continues with other methods.
      It may also, in non-standard usage, be helpful in heuristic search methods like 
      "Nelder-Mead" to avoid parameter
      regions that are unwanted or inadmissible. It is inappropriate for gradient methods.
      Default is (0.5)*.Machine$double.xmax }
  }

  There may be \code{control} elements that apply only to some of the methods. Using these
  may or may not "work" with \code{opm()}, and errors may occur with methods for which 
  the controls have no meaning. 
  However, it should be possible to call the underlying \code{optimr()} function with 
  these method-specific controls.

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Note that no other attributes of \code{par}
  are copied over. (We have not verified this as at 2009-07-29.)
}
\value{

   If there are \code{npar} parameters, then the result is a dataframe having one row
   for each method for which results are reported, using the method as the row name,
   with columns

   \code{par_1, .., par_npar, value, fevals, gevals, niter, convergence, kkt1, kkt2, xtimes}

  where
  \describe{
  \item{par_1}{ .. }
  \item{par_npar}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{fevals}{The number of calls to \code{fn}. NOT reported for method \code{lbfgs}.}
  \item{gevals}{The number of calls to \code{gr}. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.  NOT reported for method \code{lbfgs}.}
  \item{convergence}{An integer code. \code{0} indicates successful
    convergence. Various methods may or may not return sufficient information
	to allow all the codes to be specified. An incomplete list of codes includes
    \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
      \item{\code{20}}{indicates that the initial set of parameters is inadmissible, that is,
	that the function cannot be computed or returns an infinite, NULL, or NA value.}
      \item{\code{21}}{indicates that an intermediate set of parameters is inadmissible.}
      \item{\code{10}}{indicates degeneracy of the Nelder--Mead simplex.}
      \item{\code{51}}{indicates a warning from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{52}}{indicates an error from the \code{"L-BFGS-B"}
      method; see component \code{message} for further details.}
      \item{\code{9998}}{indicates that the method has been called with a NULL 'gr'
         function, and the method requires that such a function be supplied.}
      \item{\code{9999}}{indicates the method has failed.}
    }
  }
  \item{kkt1}{A logical value returned TRUE if the solution reported has a ``small'' gradient.}
  \item{kkt2}{A logical value returned TRUE if the solution reported appears to have a 
  positive-definite Hessian.}
  \item{xtimes}{The reported execution time of the calculations for the particular method.}
  }


The attribute "details" to the returned answer object contains information,
if computed, on the gradient (\code{ngatend}) and Hessian matrix (\code{nhatend}) 
at the supposed optimum, along with the eigenvalues of the Hessian (\code{hev}), 
as well as the \code{message}, if any, returned by the computation for each \code{method},
which is included for each row of the \code{details}. 
If the returned object from optimx() is \code{ans}, this is accessed 
via the construct
    \code{attr(ans, "details")}

This object is a  matrix based on a list so that if ans is the output of optimx
then attr(ans, "details")[1, ] gives the first row and 
attr(ans,"details")["Nelder-Mead", ] gives the Nelder-Mead row. There is 
one row for each method that has been successful 
or that has been forcibly saved by save.failures=TRUE. 

There are also attributes
   \describe{
   \item{maximize}{to indicate we have been maximizing the objective}
   \item{npar}{to provide the number of parameters, thereby facilitating easy
        extraction of the parameters from the results data frame}
   \item{follow.on}{to indicate that the results have been computed sequentially,
        using the order provided by the user, with the best parameters from one
        method used to start the next. There is an example (\code{ans9}) in 
        the script \code{ox.R} in the demo directory of the package.}
   }
}
\note{
  Most methods in \code{optimx} will work with one-dimensional \code{par}s, but such
  use is NOT recommended. Use \code{\link{optimize}} or other one-dimensional methods instead.

  There are a series of demos available. Once the package is loaded (via \code{require(optimx)} or
  \code{library(optimx)}, you may see available demos via 

  demo(package="optimx")

  The demo 'brown_test' may be run with the command
  demo(brown_test, package="optimx")

  The package source contains several functions that are not exported in the
  NAMESPACE. These are 
  \describe{
  \item{\code{optimx.setup()}}{ which establishes the controls for a given run;}
  \item{\code{optimx.check()}}{ which performs bounds and gradient checks on
      the supplied parameters and functions;}
  \item{\code{optimx.run()}}{which actually performs the optimization and post-solution
      computations;}
  \item{\code{scalechk()}}{ which actually carries out a check on the relative scaling
      of the input parameters.}
  }

Knowledgeable users may take advantage of these functions if they are carrying
out production calculations where the setup and checks could be run once.

}
\source{

See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

}
\references{

 See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

 Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

 Nash JC (2014). On Best Practice Optimization Methods in R., 
        \emph{Journal of Statistical Software}, 60(2), 1-14.,
        URL http://www.jstatsoft.org/v60/i02/.

}

\seealso{
  \code{\link[BB]{spg}}, \code{\link{nlm}}, \code{\link{nlminb}},
  \code{\link[minqa]{bobyqa}}, 
  \code{\link[ucminf]{ucminf}}, 
  \code{\link[dfoptim]{nmkb}},
  \code{\link[dfoptim]{hjkb}}.
  \code{\link{optimize}} for one-dimensional minimization;
  \code{\link{constrOptim}} or \code{\link[BB]{spg}} for linearly constrained optimization.
}

\examples{
require(graphics)
cat("Note possible demo(ox) for extended examples\n")


## Show multiple outputs of optimx using all.methods
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
## One generalization of the Rosenbrock banana valley function (n parameters)
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
        return(fval)
}

genrose.g <- function(x, gs=NULL){
# vectorized gradient for genrose.f
# Ravi Varadhan 2009-04-03
	n <- length(x)
        if(is.null(gs)) { gs=100.0 }
	gg <- as.vector(rep(0, n))
	tn <- 2:n
	tn1 <- tn - 1
	z1 <- x[tn] - x[tn1]^2
	z2 <- 1 - x[tn]
	gg[tn] <- 2 * (gs * z1 - z2)
	gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
	return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
   if(is.null(gs)) { gs=100.0 }
	n <- length(x)
	hh<-matrix(rep(0, n*n),n,n)
	for (i in 2:n) {
		z1<-x[i]-x[i-1]*x[i-1]
		z2<-1.0-x[i]
                hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
                hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
                hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
                hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
	}
        return(hh)
}

startx<-4*seq(1:10)/3.
ans8<-opm(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h,
   method="ALL", control=list(save.failures=TRUE, trace=0), gs=10)
# Set trace=1 for output of individual solvers
ans8
ans8[, "gevals"]
ans8["spg", ]
summary(ans8, par.select = 1:3)
summary(ans8, order = value)[1, ] # show best value
head(summary(ans8, order = value)) # best few
## head(summary(ans8, order = "value")) # best few -- alternative syntax

## order by value.  Within those values the same to 3 decimals order by fevals.
## summary(ans8, order = list(round(value, 3), fevals), par.select = FALSE)
summary(ans8, order = "list(round(value, 3), fevals)", par.select = FALSE)

## summary(ans8, order = rownames, par.select = FALSE) # order by method name
summary(ans8, order = "rownames", par.select = FALSE) # same

summary(ans8, order = NULL, par.select = FALSE) # use input order
## summary(ans8, par.select = FALSE) # same

}
\keyword{nonlinear}
\keyword{optimize}
