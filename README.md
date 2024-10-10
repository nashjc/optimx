# optimx

optimx is an R package that extends and enhances the optim() function of base R,
in particular by unifying the call to many solvers.
This Gitlab project has been established to document optimx and its development.
A more detailed history of the package is given at the bottom of this document.

## Installation

From within R, the command 

    install.packages("optimx")

will install the optimx from the CRAN repository https://cran.r-project.org/

of which a copy is here on github.com/nashjc/. Thus if the devtools package is
installed in R on your computer, you can issue the commands

    library("devtools")
    install_github("nashjc/optimx")

The test version of optimx on R-forge in the optimizer project 
(https://r-forge.r-project.org/R/?group_id=395) is now deprecated.

## Main functions in optimx

### optimr()

optimr() allows any of the solvers linked to the package to be called via the call

  result <- optimr(par, fn, gr, hess, method="a_solver", lower=lo, upper=up)
          
where
  - par is a vector of the initial parameter values for the function
  - fn is the objective function to be minimized
  - gr is a function to compute the gradient of fn (or NULL if none is provided)
  - hess is a function to compute the hessian of fn (or NULL if gr=NULL or none is provided)
  - method provides the name of the solver, here given a token name
  - lower and upper provide bounds on the parameters for solvers that can use them

See the man page for other arguments that allow various controls or extra features.

### opm()

opm() allows multiple solvers to
be called via a single statement so that their performance may be compared.
It is INEFFICIENT, WASTEFUL, and strongly DISCOURAGED to use opm() as a 
general-purpose optimization tool. The call is similar to that of optimr()
but method is now a vector of names of solvers. See the man page for more
details.

Usage:

result<-opm(par, fn, gr, hess, lower=lo, upper=up, method=c("solver1", "solver2", "solver3"))

with arguments as for optimr() except for method, which is now a vector of character names.

## Solvers provided within the optimx package

The package contains several solver functions. See specific man pages for details. Some of these
are intended to be called via optimr() and NOT independently.

 - ncg() is a conjugate gradient minimizer based on the Dai/Yuan approach. 
 - Rcgmin() is an earlier implementation of ncg() which may give slightly different results and
   is provided to maintain backward compatibility. 
 - Rcgminb() is the bounds-constrained code called by Rcgmin()
 - Rcgminu() is the unconstrained  code called by Rcgmin()
 - nvm() is a Fletcher variable metric minimizer similar to the 'BFGS' method of optim()
   (also available via optimr())
 - Rvmmin() is an earlier implementation of nvm() which may give slightly different results and
   is provided to maintain backward compatibility.
 - Rvmminb() is the bounds-constrained code called by Rvmmin()
 - Rvmminu() is the unconstrained  code called by Rvmmin()
 - snewtm() is a didactic safeguarded Newton method with a Marquardt stabilization of the Hessian
 - snewton() is a didactic safeguarded Newton method with a backtracking line search
 - tn() is an unconstrained truncated Newton method of Stephen Nash
 - tnbc() is a bounds-constrained truncated Newton method of Stephen Nash
 - Rtnmin-package contains functions needed by tn() or tnbc()
 - hjn() is a didactic implementation of the Hooke and Jeeves pattern search

## Solver characteristics

The table includes methods from other packages or base-R as well as solvers in optimx.

| solver        | source    | gradient | hessian  | bounds| algorithm  class |  Notes                 |
| ------------- | --------- |  ------- |  ------- | ----- | ---------------- | ---------------------- |
| anms          | pracma    |          |          |       | polytope         |                        |
| BFGS          | optim()   | optional |          |       | variable metric  |                        |
| bobyqa        | minqa     |          |          |  yes  | model & descend  |                        |
| CG            | optim()   | optional |          |       | conjugate grad   | Deprecated             |
| hjkb          | dfoptim   |          |          |  yes  | Hooke&Jeeves     |                        |
| hjn           | optimx    |          |          |  yes  | Hooke&Jeeves     |                        |
| lbfgs         | lbfgs     | optional |          |       | LM quasi Newton  |                        |
| L-BFGS-B      | optim()   | optional |          |  yes  | LM quasi Newton  |                        |
| lbfgsb3c      | lbfgsb3c  | optional |          |  yes  | LM quasi Newton  |                        |
| mla           | marqLevAlg| optional | optional |       | Newton-like      |                        |
| ncg           | optimx    | required |          |  yes  | conjugate grad   |                        |
| Nelder-Mead   | optim()   |          |          |       | polytope         |                        |
| newuoa        | minqa     |          |          |       | model & descend  |                        |
| nlm           | base-R    | optional | optional |       | Newton-like      |                        |
| nlminb        | base-R    | optional | optional |  yes  | multiple         |                        |
| nlnm          | nloptr    |          |          |  yes  | polytope         |                        |
| nmkb          | dfoptim   |          |          |  yes  | polytope         |                        |
| nvm           | optimx    | required |          |  yes  | variable metric  |                        |
| pracmanm      | pracma    |          |          |       | polytope         |                        |
| Rcgmin        | optimx    | required |          |  yes  | conjugate grad   | legacy compatibility   |
| Rtnmin        | optimx    | required |          |  yes  | trunc. Newton    |                        |
| Rvmmin        | optimx    | required |          |  yes  | variable metric  | legacy compatibility   |
| slsqp         | nloptr    | optional |          |  yes  | quad programming |                        |
| snewtm        | optimx    | required | required |       | Newton-like      | Didactic               |
| snewton       | optimx    | required | required |       | Newton-like      | Didactic               |
| spg           | BB        | optional |          |  yes  | proj. gradient   |                        |
| subplex       | subplex   |          |          |       | polytope         |                        |
| tnewt         | nloptr    | optional |          |  yes  | trunc. Newton    |                        |
| ucminf        | ucminf    | optional |          |       | variable metric  |                        |
| uobyqa        | minqa     |          |          |       | model & descend  | Use newuoa             |
|---------------|-----------|----------|----------|-------|------------------|------------------------|

Note that there are other versions of some of these algorithms
in different packages, e.g., package nloptr has a bobyqa() function.

## Adding solvers

It is relatively easy to modify optimr.R and ctrldefault.R to add a new solver. See the vignette
"Using and extending the optimx package" (file Extend-optimx.pdf).

## Additional functions within optimx

The following files provide tools to support the use of the optimx package
and other optimization tasks. See the man pages for details on usage.

- axsearch.R: Perform axial search around a supposed minimum and provide diagnostics
- bmchk.R: Check bounds and masks for parameter constraints
- bmstep.R: Compute the maximum step along a search direction in presence of bounds
- checksolver.R: Test if requested solver is present
- ctrldefault.R: Set control defaults
- fnchk.R: Run tests, where possible, on user objective function
- gHgen.R: Generate gradient and Hessian for a function at given parameters
- gHgenb.R: Generate gradient and Hessian for a function at given parameters (with bounds constraints)
- grchk.R: Run tests, where possible, on user objective function and (optionally) gradient and hessian
- grback.R: Backward difference numerical gradient approximation
- grcentral.R: Central difference numerical gradient approximation 
- grfwd.R: Forward difference numerical gradient approximation
- grnd.R: A reorganization of the call to numDeriv grad() function
- grpracma.R: A reorganization of the call to the pracma grad() function
- hesschk.R: Run tests, where possible, on user objective function and (optionally) gradient and hessian
- kktchk.R: Check Kuhn Karush Tucker conditions for a supposed function minimum
- multistart.R: Allows for multiple sets of optimization starting parameters to be tried in a single call
- opm2optimr.R: Extract optim() solution for one method of opm() result
- optimx.R: The original optimx() function from 2010 for backward compatibility. This calls optimx.setup(),
  optimx.check() and optimx.run().
- polyopt.R: Allows solvers to be applied sequentially, with control limits on computational effort for
  each solver
- proptimr.R: Compact display of an optimr() result object 
- scalechk.R: Check scale of initial parameters and bounds supplied for an optimization problem 
- optchk.R: Attempts to check user-supplied objective function (optionally also gradient and hessian)
- optimx-package.R: Provides local versions of summary() and coef() functions and some other code.
- zzz.R: Startup actions. At time of writing, used for tn() and tnbc().

## Fixed parameters via equal bounds (masks)

Some solvers, in particular ncg() and nvm() can explicitly deal with
lower and upper bounds constraints that are set equal for some parameters. Such fixed
parameters can be useful when a parameter is not generally altered
to try to optimize the function, but may be allowed into the optimization
in the future. Such constraints are sometimes called "masks". The initial
parameter value must be consistent with the bounds settings.

Use of equal lower and upper bounds may be possible for other solvers that
handle bounds constraints, but users are urged to test this possibility before
trusting results.


## History of optimx

The optimx package was first placed on CRAN in 2009 by John Nash and Ravi
Varadhan. The purpose was to allow a single syntax to call one or several 
different optimization solvers with a single optimization problem. The 
package was quite successful, and the optimx() function within the package
is still present in the current package with some minor corrections. However,
the structure of the optimx() function made it difficult to extend and
maintain. Moreover, the output was slightly altered from that of the base-R
optim() function. 

The optimr and optimrx (experimental) packages were an attempt to resolve
the maintenance issue and to address some deficiencies in optimx(). In 
particular, the optimr() function was introduced to mirror the base-R
optim() function but with more solvers. The multiple-solver feature of
the optimx() function was replaced with opm(), while polyalgorithm and
multiple initial parameter sets (starts) were handled with polyopt() and
multistart() respectively, avoiding some complicated option possibilities
in optimx(). Moreover, the optimr() function now permits all solvers to 
use parameter scaling via the control element (actually a vector) parscale.

In 2016/2017 it became clear that there were a number of optimization tools 
that could be consolidated into a single packages. The current optimx package
replaces:

   - optimx
   - optimr
   - Rvmmin
   - Rcgmin
   - Rtnmin
   - optextras
   
It also includes two new Newton-like solvers, a didactic Hooke and Jeeves
solver, and a simple interface to various approximate gradient routines.
As and when it is sensible to do so, the subsumed packages will be archived.
   
Note that some development versions of packages remain on R-forge. Use of these 
is at your own risk. 


### Updated 2024-10-03
