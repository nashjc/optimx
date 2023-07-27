\name{optimx-package}
\alias{optimx-package}
\docType{package}
\title{
A replacement and extension of the optim() function, plus various 
optimization tools
}
\description{
\code{\link{optimx}}
provides a replacement and extension of the \code{link{optim()}} function to 
unify and streamline optimization capabilities in R for smooth, possibly box
constrained functions of several or many parameters

The three functions ufn, ugr and uhess wrap corresponding user functions fn, gr, and 
hess so that these functions can be executed safely (via try()) and also so parameter or
function scaling can be applied. The wrapper
functions also allow for maximization of functions (via minimization of the negative of
the function) using the logical parameter \code{maximize}.

There are three test functions, fnchk, grchk, and hesschk, to allow the user 
function to be tested for validity and correctness. However, no set of tests is 
exhaustive, and extensions and improvements are welcome. The package 
\code{numDeriv} is used for generation of numerical approximations to 
derivatives.

}
\details{

Index:
\preformatted{
axsearch            Perform an axial search optimality check
bmchk               Check bounds and masks for parameter constraints
bmstep              Compute the maximum step along a search direction.
checksolver         Checks if method is available in allmeth
ctrldefault         Sets the default values of elements of the control() list
fnchk               Test validity of user function
gHgen               Compute gradient and Hessian as a given 
                    set of parameters
gHgenb              Compute gradient and Hessian as a given 
                    set of parameters appying bounds and masks
grback              Backward numerical gradient approximation
grcentral           Central numerical gradient approximation
grchk               Check that gradient function evaluation 
                    matches numerical gradient
grfwd               Forward numerical gradient approximation
grnd                Gradient approximation using \code{numDeriv}
hesschk             Check that Hessian function evaluation 
                    matches numerical approximation
hjn                 A didactic example code of the Hooke and Jeeves algorithm
kktchk              Check the Karush-Kuhn-Tucker optimality conditions
scalechk            Check scale of initial parameters and bounds
multistart          Try a single method with multiple starting parameter sets
opm                 Wrapper that allows multiple minimizers to be applied to a
                    given objective function 
optimr              Wrapper that allows different (single) minimizers to be 
                    applied to a given objective function using a common syntax
                    like that of optim()
optimx              Wrapper that allows multiple minimizers to be applied to a
                    given objective function. Complexity of the code maked this
                    function difficult to maintain, and opm() is the suggested
                    replacement, but optimx() is retained for backward 
                    compatibility.
optimx.check        a component of optimx()
optimx-package      a component of optimx()
optimx.run          a component of optimx()
optimx.setup        a component of optimx()
optsp               An environment to hold some globally useful items
                    used by optimization programs. Created on loading package
                    with zzz.R
polyopt             Allows sequential application of methods to a given problem.
proptimr            compact output of optimr() result object             
Rcgminb             Bounds constrained conjugate gradients minimization
Rcgmin              Conjugate gradients minimization
Rcgminu             Unconstrained conjugate gradients minimization
Rtnmin-package      Internal functions for the S.G.Nash truncated newton method 
Rvmminb             Bounds constrained variable metric minimization method
Rvmmin              Variable metric minimization method
Rvmminu             Unconstrained variable metric minimization method
snewtonm            Demonstration Newton-Marquardt minimization method
snewton             Demonstration safeguarded Newton minimization method
snewtonmb           Bounds constrained safeguarded Newton method
tnbc                Bounds constrained truncated Newton method
tn                  Unconstrained truncated Newton method
}
}
\author{
John C Nash <nashjc@uottawa.ca> and Ravi Varadhan <RVaradhan@jhmi.edu>

Maintainer: John C Nash <nashjc@uottawa.ca>
}
\references{
  Nash, John C. and Varadhan, Ravi (2011) Unifying Optimization Algorithms 
          to Aid Software System Users: optimx for R, Journal of Statistical
	        Software, publication pending.
}
\keyword{ package }
\keyword{ optimization }

