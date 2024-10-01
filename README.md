# optimx

optimx is an R package that extends and enhances the optim() function of base R. 
This Gitlab project has been established to document optimx and its development.
A more detailed history of the package is given after the links to it.

## Purpose

There are many function minimization and optimization packages associated with R.
Generally, "optimization" includes constraints, but in the current context we
use it to mean function minimization or maximization with at most bounds 
constraints. We discuss all problems as minimizations, though the tools generally
offer a setting to maximize a function.
Unfortunately, each may have particular syntax and other requirements. Thus the
primary purpose of optimx, and in particular the optimr() function within the 
package, is intended to allow a common calling syntax. In large measure, this 
syntax mirrors that of the base-R optim() function. For that reason, some features 
of the solvers callable by optimx may not be available with optimr().

A second goal of optimx is, via the function opm(), to allow multiple solvers to
be called via a single statement so that their performance may be compared. When
a user has a task that will be repeated frequently, this permits a reasonable 
choice of efficient solver for the problem at hand. Some users have thought to
use opm() generally with all, or many, of the available solvers, choosing the
"best" solution after the fact. This is INEFFICIENT, WASTEFUL, and strongly
DISCOURAGED.

A third role of optimx is to provide several solvers and some infrastructure
tools for function minimization. These were previously available in packages
Rvmmin, Rcgmin, Rtnmin and optextras. However, the overhead of updating each
of these packages separately became overly burdensom for the maintainer, and
now all these features are included within optimx. We urge developers of other
packages that called the separate tools to adjust their code to call optimx. 
This is largely a matter of simply adjusting require() or library() entries 
or similar simple changes.


## Where to find optimx

The official optimx package is on CRAN. 

https://cran.r-project.org/

A test version has been on R-forge in the optimizer package for some time, 
but there have been concerns that R-forge is not up to date:

https://r-forge.r-project.org/R/?group_id=395

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

## Categorization of solvers

It is helpful to have a general overview of the solvers available in optimx.

Possibly the most important consideration is whether the solvers require or
assume that gradients are available, or whether we may describe the solver
as a direct search. Note that numerical derivatives do allow gradient or 
Newton-like methods to be applied without explicit code for the gradient
or Hessian. However, this is often inefficient, though timings can be useful
to help choose appropriate solvers.

The solvers that do not require gradients are named "Nelder-Mead", "newuoa", 
"bobyqa", "uobyqa", "nmkb", "hjkb", "hjn", "subplex", "anms", "pracmanm", and "nlnm".
Several of these are variants on the Nelder-Mead polytope (or simplex) search. 
The Powell methods that use an "approximate and descend" approach are "newuoa", 
"bobyqa", and "uobyqa". However, it has recently been reported that these may
have bugs. See T. M. Ragonneau and Z. Zhang, PDFO: a cross-platform package for 
Powell's derivative-free optimization solvers, arXiv:2302.13246, 2023

Of methods that require gradients, those using ideas from Newton's method 
require storage of some approximation to the Hessian. Methods that can use
an explicit Hessian are "nlm", "nlminb", "snewton", "snewtonm", and "snewtm".
The last three require the Hessian and are mainly didactic or experimental.
The other two can use an explicit Hessian but do not require it. Timings are
recommended to decide on usage details. 

Clever ideas to approximate the Hessian without explicit calculation of it may
still require storage for the matrix. When the problem dimension is large, there
may be memory implications, though modern computers generally have plenty of 
fast memory. "BFGS", "nlm", "nlminb", "Rvmmin", "ucminf", "nvm", "mla", and
"slsqp" use explicit storage of an approximate Hessian.

Avoiding the explicit matrix for the approximate Hessian leads to methods that
are called "limited memory" or "conjugate gradient" or "truncated Newton" methods.
"CG", "L-BFGS-B", "lbfgsb3c", "Rcgmin", "Rtnmin", "spg", "lbfgs", "ncg", and "tnewt"
fall in this category.

A separate categorization describes which methods allow bounds constraints.
"L-BFGS-B", "nlminb", "lbfgsb3c", "Rcgmin", "Rtnmin", "nvm", 
"Rvmmin", "bobyqa", "nmkb", "hjkb", "hjn", "snewtonm", "ncg", "slsqp", "tnewt", "nlnm", 
"snewtm", and "spg" provide for such constraints. Some, but not all, of these
will allow upper and lower bounds on some parameters to be equal, thereby fixing a parameter
or establishing a "mask". This can be very useful for situations where a parameter has
a generally accepted value and we choose not to try to optimizat over it. Later the
constraint can be relaxed. In particular "ncg" and "nvm" have been set up to explicitly
handle masks.

Note that some solvers are updated versions of existing codes. "nvm" updates "Rvmmin" and "ncg"
updates "Rcgmin", though there are some situations where the older codes appear to perform
better. Moreover, users can get upset if tests do not reproduce earlier results, so the older
solvers have been left available. 

## Adding solvers

It is relatively easy to modify optimr.R and ctrldefault.R to add a new solver. See the vignette
"Using and extending the optimx package" (file Extend-optimx.pdf).

### Updated 2024-10-01
