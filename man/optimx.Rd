% File src/library/stats/man/optimx.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{optimx}
\alias{optimx}
\alias{[.optimx}
\alias{as.data.frame.optimx}
\encoding{UTF-8}
\title{General-purpose optimization}
\concept{minimization}
\concept{maximization}
\description{
  General-purpose optimization wrapper function that calls other
  R tools for optimization, including the existing optim() function.
  \code{optimx} also tries to unify the calling sequence to allow
  a number of tools to use the same front-end. These include 
  \code{spg} from the BB package, \code{ucminf}, \code{nlm}, and 
  \code{nlminb}. Note that 
  optim() itself allows Nelder--Mead, quasi-Newton and 
  conjugate-gradient algorithms as well as box-constrained optimization 
  via L-BFGS-B. Because SANN does not return a meaningful convergence code
  (conv), optimx() does not call the SANN method.
  
  Note that package \code{optimr} allows solvers to be called individually
  by the \code{optim()} syntax, with the \code{parscale} 
  control to scale parameters applicable to all methods. However, 
  running multiple methods, or using the \code{follow.on} capability
  has been moved to separate routines in the \code{optimr} package.
  
  Cautions: 
  
  1) Using some control list options with different or multiple methods 
  may give unexpected results. 
  
  2) Testing the KKT conditions can take much longer than solving the
  optimization problem, especially when the number of parameters is large
  and/or analytic gradients are not available. Note that the default for
  the control \code{kkt} is TRUE.
}

\usage{
optimx(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf, 
            method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
            control=list(),
             ...)
}
\arguments{
 \item{par}{a vector of initial values for the parameters 
   for which optimal values are to be found. Names on the elements
   of this vector are preserved and used in the results data frame.}  
 \item{fn}{A function to be minimized (or maximized), with first
   argument the vector of parameters over which minimization is to take
   place.  It should return a scalar result.}
 \item{gr}{A function to return (as a vector) the gradient for those methods that 
   can use this information.

   If 'gr' is \code{NULL}, a finite-difference approximation will be used.
   An open question concerns whether the SAME approximation code used for all methods, 
   or whether there are differences that could/should be examined? }

 \item{hess}{A function to return (as a symmetric matrix) the Hessian of the objective 
   function for those methods that can use this information.}
 \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can
   handle box (or bounds) constraints.}
 \item{method}{A list of the methods to be used. 
       Note that this is an important change from optim() that allows
       just one method to be specified. See \sQuote{Details}.}
 \item{itnmax}{If provided as a vector of the same length as the list of methods \code{method}, 
	gives the maximum number of iterations or function values for the corresponding 
	method. If a single number is provided, this will be used for all methods. Note that
	there may be control list elements with similar functions, but this should be the
	preferred approach when using \code{optimx}.}
 \item{hessian}{A logical control that if TRUE forces the computation of an approximation 
       to the Hessian at the final set of parameters. If FALSE (default), the hessian is
       calculated if needed to provide the KKT optimality tests (see \code{kkt} in
       \sQuote{Details} for the \code{control} list).
       This setting is provided primarily for compatibility with optim().}
 \item{control}{A list of control parameters. See \sQuote{Details}.}
 \item{\dots}{For \code{optimx} further arguments to be passed to \code{fn} 
    and \code{gr}; otherwise, further arguments are not used.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$maximize} is TRUE. The original optim() function allows
  \code{control$fnscale} to be set negative to accomplish this. DO NOT
  use both methods. 

  Possible method codes at the time of writing are 'Nelder-Mead', 'BFGS',
  'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa',
  'nmkb', 'hjkb', 'Rcgmin', or 'Rvmmin'.

  The default methods for unconstrained problems (no \code{lower} or
  \code{upper} specified) are an implementation of the Nelder and Mead
  (1965) and a Variable Metric method based on the ideas of Fletcher
  (1970) as modified by him in conversation with Nash (1979). Nelder-Mead
  uses only function values and is robust but relatively slow.  It will 
  work reasonably well for non-differentiable functions. The Variable
  Metric method, \code{"BFGS"} updates an approximation to the inverse
  Hessian using the BFGS update formulas, along with an acceptable point
  line search strategy. This method appears to work best with analytic
  gradients. (\code{"Rvmmmin"} provides a box-constrained version of this
  algorithm.

  If no \code{method} is given, and there are bounds constraints provided,
  the method is set to \code{"L-BFGS-B"}.

  Method \code{"CG"} is a conjugate gradients method based on that by
  Fletcher and Reeves (1964) (but with the option of Polak--Ribiere or
  Beale--Sorenson updates). The particular implementation is now dated,
  and improved yet simpler codes are being implemented (as at June 2009),
  and furthermore a version with box constraints is being tested.
  Conjugate gradient methods will generally be more fragile than the 
  BFGS method, but as they do not store a matrix they may be successful 
  in much larger optimization problems.

  Method \code{"L-BFGS-B"} is that of Byrd \emph{et. al.} (1995) which
  allows \emph{box constraints}, that is each variable can be given a lower
  and/or upper bound. The initial value must satisfy the constraints.
  This uses a limited-memory modification of the BFGS quasi-Newton
  method. If non-trivial bounds are supplied, this method will be
  selected, with a warning.

  Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.

  Function \code{fn} can return \code{NA} or \code{Inf} if the function
  cannot be evaluated at the supplied value, but the initial value must
  have a computable finite value of \code{fn}. However, some methods, of
  which \code{"L-BFGS-B"} is known to be a case, require that the values
  returned should always be finite.

  While \code{optim} can be used recursively, and for a single parameter
  as well as many, this may not be true for \code{optimx}. \code{optim}
  also accepts a zero-length \code{par}, and just evaluates the function 
  with that argument.

  Method \code{"nlm"} is from the package of the same name that implements
  ideas of Dennis and Schnabel (1983) and Schnabel et al. (1985). See nlm()
  for more details.

  Method \code{"nlminb"} is the package of the same name that uses the
  minimization tools of the PORT library.  The PORT documentation is at 
  <URL: http://netlib.bell-labs.com/cm/cs/cstr/153.pdf>. See nlminb()
  for details. (Though there is very little information about the methods.)

  Method \code{"spg"} is from package BB implementing a spectral projected 
  gradient method for large-scale optimization with simple constraints due
  R adaptation, with significant modifications, by Ravi Varadhan,
  Johns Hopkins University (Varadhan and Gilbert, 2009), from the original
  FORTRAN code of Birgin, Martinez, and Raydan (2001). 

  Method \code{"Rcgmin"} is from the package of that name. It implements a
  conjugate gradient algorithm with the Dai and Yuan (2001) update and also 
  allows bounds constraints on the parameters. (Rcgmin also allows mask 
  constraints -- fixing individual parameters.) 

  Methods \code{"bobyqa"}, \code{"uobyqa"} and \code{"newuoa"} are from the 
  package \code{"minqa"} which implement optimization by quadratic approximation
  routines of the similar names due to M J D Powell (2009). See package minqa 
  for details. Note that \code{"uobyqa"} and \code{"newuoa"} are for 
  unconstrained minimization, while \code{"bobyqa"} is for box constrained
  problems. While \code{"uobyqa"} may be specified, it is NOT part of the 
  \code{all.methods = TRUE} set.

  The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{trace}}{Non-negative integer. If positive,
      tracing information on the
      progress of the optimization is produced. Higher values may
      produce more tracing information: for method \code{"L-BFGS-B"}
      there are six levels of tracing. trace = 0 gives no output 
      (To understand exactly what these do see the source code: higher 
      levels give more detail.)}
    \item{\code{follow.on }}{ = TRUE or FALSE. If TRUE, and there are multiple 
      methods, then the last set of 
      parameters from one method is used as the starting set for the next. }
    \item{\code{save.failures}}{ = TRUE if we wish to keep "answers" from runs 
      where the method does not return convcode==0. FALSE otherwise (default).}
    \item{\code{maximize}}{ = TRUE if we want to maximize rather than minimize 
      a function. (Default FALSE). Methods nlm, nlminb, ucminf cannot maximize a
      function, so the user must explicitly minimize and carry out the adjustment
      externally. However, there is a check to avoid
      usage of these codes when maximize is TRUE. See \code{fnscale} below for 
      the method used in \code{optim} that we deprecate.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkt}}{=FALSE if we do NOT want to test the Kuhn, Karush, Tucker
      optimality conditions. The default is TRUE. However, because the Hessian
      computation may be very slow, we set \code{kkt} to be FALSE if there are 
      more than than 50 parameters when the gradient function \code{gr} is not 
      provided, and more than 500
      parameters when such a function is specified. We return logical values \code{KKT1}
      and \code{KKT2} TRUE if first and second order conditions are satisfied approximately.
      Note, however, that the tests are sensitive to scaling, and users may need
      to perform additional verification. If \code{kkt} is FALSE but \code{hessian}
      is TRUE, then \code{KKT1} is generated, but \code{KKT2} is not.}
    \item{\code{all.methods}}{= TRUE if we want to use all available (and suitable)
      methods.}
    \item{\code{kkttol}}{= value to use to check for small gradient and negative
      Hessian eigenvalues. Default = .Machine$double.eps^(1/3) }
    \item{\code{kkt2tol}}{= Tolerance for eigenvalue ratio in KKT test of positive 
      definite Hessian. Default same as for kkttol }
    \item{\code{starttests}}{= TRUE if we want to run tests of the function and 
      parameters: feasibility relative to bounds, analytic vs numerical gradient, 
      scaling tests, before we try optimization methods. Default is TRUE.}
    \item{\code{dowarn}}{= TRUE if we want warnings generated by optimx. Default is 
      TRUE.}
    \item{\code{badval}}{= The value to set for the function value when try(fn()) fails.
      Default is (0.5)*.Machine$double.xmax }
    \item{\code{usenumDeriv}}{= TRUE if the \code{numDeriv} function \code{grad()} is
      to be used to compute gradients when the argument \code{gr} is NULL or not supplied.}
  }

  The following \code{control} elements apply only to some of the methods. The list
  may be incomplete. See individual packages for details. 

  \describe{

    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} and \code{gr} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}. For methods from the set in
      \code{optim()}. Note potential conflicts with the control \code{maximize}.}
    \item{\code{parscale}}{A vector of scaling values for the parameters.
	Optimization is performed on \code{par/parscale} and these should be
	comparable in the sense that a unit change in any element produces
	about a unit change in the scaled value.For \code{optim}.}
    \item{\code{ndeps}}{A vector of step sizes for the finite-difference
      approximation to the gradient, on \code{par/parscale}
      scale. Defaults to \code{1e-3}. For \code{optim}.}
    \item{\code{maxit}}{The maximum number of iterations. Defaults to
      \code{100} for the derivative-based methods, and
      \code{500} for \code{"Nelder-Mead"}.}
    \item{\code{abstol}}{The absolute convergence tolerance. Only
      useful for non-negative functions, as a tolerance for reaching zero.}
    \item{\code{reltol}}{Relative convergence tolerance.  The algorithm
      stops if it is unable to reduce the value by a factor of
      \code{reltol * (abs(val) + reltol)} at a step.  Defaults to
      \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}. For \code{optim}.}
    \item{\code{alpha}, \code{beta}, \code{gamma}}{Scaling parameters
      for the \code{"Nelder-Mead"} method. \code{alpha} is the reflection
      factor (default 1.0), \code{beta} the contraction factor (0.5) and
      \code{gamma} the expansion factor (2.0).}
    \item{\code{REPORT}}{The frequency of reports for the \code{"BFGS"} and
      \code{"L-BFGS-B"} methods if \code{control$trace}
      is positive. Defaults to every 10 iterations for \code{"BFGS"} and
      \code{"L-BFGS-B"}.}
    \item{\code{type}}{for the conjugate-gradients method. Takes value
      \code{1} for the Fletcher--Reeves update, \code{2} for
      Polak--Ribiere and \code{3} for Beale--Sorenson.}
    \item{\code{lmm}}{is an integer giving the number of BFGS updates
      retained in the \code{"L-BFGS-B"} method, It defaults to \code{5}.}
    \item{\code{factr}}{controls the convergence of the \code{"L-BFGS-B"}
      method. Convergence occurs when the reduction in the objective is
      within this factor of the machine tolerance. Default is \code{1e7},
      that is a tolerance of about \code{1e-8}.}
    \item{\code{pgtol}}{helps control the convergence of the \code{"L-BFGS-B"}
      method. It is a tolerance on the projected gradient in the current
      search direction. This defaults to zero, when the check is
      suppressed.}
  }

  Any names given to \code{par} will be copied to the vectors passed to
  \code{fn} and \code{gr}.  Note that no other attributes of \code{par}
  are copied over. (We have not verified this as at 2009-07-29.)

  There are \code{[.optimx}, \code{as.data.frame.optimx}, \code{\link{coef.optimx}} 
  and \code{\link{summary.optimx}} methods available.
  
  Note: Package \code{optimr} is a derivative of this package. It was developed
  initially to overcome maintenance difficulties with the current package 
  related to avoiding confusion if some multiple options were specified together,
  and to allow the \code{optim()} function syntax to be used consistently, 
  including the \code{parscale} control. However, this package does perform
  well, and is called by a number of popular other packages.

}
\value{

   If there are \code{npar} parameters, then the result is a dataframe having one row
   for each method for which results are reported, using the method as the row name,
   with columns

   \code{par_1, .., par_npar, value, fevals, gevals, niter, convcode, kkt1, kkt2, xtimes}

  where
  \describe{
  \item{par_1}{ .. }
  \item{par_npar}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{fevals}{The number of calls to \code{fn}.}
  \item{gevals}{The number of calls to \code{gr}. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to \code{fn} to
    compute a finite-difference approximation to the gradient.}
  \item{niter}{For those methods where it is reported, the number of ``iterations''. See
    the documentation or code for particular methods for the meaning of such counts.}
  \item{convcode}{An integer code. \code{0} indicates successful
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
  \item{\code{scalecheck()}}{ which actually carries out a check on the relative scaling
      of the input parameters.}
  }

Knowledgeable users may take advantage of these functions if they are carrying
out production calculations where the setup and checks could be run once.

}
\source{

See the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

}
\references{

 See also  the manual pages for \code{optim()} and the packages the DESCRIPTION \code{suggests}.

  Byrd RH, Lu P, Nocedal J (1995) A Limited Memory Algorithm for Bound Constrained Optimization,
      SIAM Journal on Scientific Computing, 16 (5), 1190--1208.

  Y. H. Dai and Y. Yuan, (2001) An Efficient Hybrid Conjugate Gradient Method for Unconstrained
  Optimization, Annals of Operations Research, 103, pp 33--47, 
  URL http://dx.doi.org/10.1023/A:1012930416777.

  Dennis JE and Schnabel RB (1983) Numerical Methods for Unconstrained Optimization and Nonlinear Equations,
     Englewood Cliffs NJ: Prentice-Hall.


   Fletcher R (1970) A New Approach to Variable Metric Algorithms, Computer Journal,
     13 (3), 317-322.

  Nash JC, and Varadhan R (2011). Unifying Optimization Algorithms to Aid Software System Users: 
    \bold{optimx} for R., \emph{Journal of Statistical Software}, 43(9), 1-14.,  
     URL http://www.jstatsoft.org/v43/i09/.

  Nash JC (2014). On Best Practice Optimization Methods in R., 
        \emph{Journal of Statistical Software}, 60(2), 1-14.,
        URL http://www.jstatsoft.org/v60/i02/.
        
  Nelder JA and Mead R (1965) A Simplex Method for Function Minimization,
     Computer Journal, 7 (4), 308--313.

  Powell MJD (2009) The BOBYQA algorithm for bound constrained optimization 
      without derivatives, http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf      

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
cat("Note demo(ox) for extended examples\n")


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
ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
   control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
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
