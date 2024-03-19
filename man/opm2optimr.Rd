\name{opm2optimr}
\alias{opm2optimr}
\encoding{UTF-8}
\title{Extract optim() solution for one method of opm() result}
\concept{minimization}
\concept{maximization}
\description{
  A function that attempts to extract a solution found for one of the
  multiple solvers used to obtain results of minimization or 
  maximinzation using the opm() multi-solver wrapper.
}
\usage{
opm2optimr(opmobj, rid)
}
\arguments{
 \item{opmobj}{the object returned by opm() from trying to solve an
  optimization problem with multiple solvers}
 \item{rid}{The identifier of the solver for which the solution should be
  extracted. This may be a "row ID" in the form of an integer or else a
  character string for the corresponding solver.}
}
\details{
  We do not extract the true "message" for a method. 
}
\value{
  A list of the following items:
 \describe{
    \item{par}{The best set of parameters found.}
    \item{value}{The value of fn corresponding to par.}
    \item{counts}{A two-element integer vector giving the number of calls to fn and gr 
            respectively. This excludes those calls needed to compute the Hessian
            even though the opm() result will have these counts}
    \item{convergence}{An integer code. 0 indicates successful completion}
    \item{message}{A character string which for optim() or optimr() may give additional 
            information returned by the optimizer, or NULL.
            Here will be "Result of conversion from opm() result"}
 }
} 
\examples{

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}
mset <- c("ncg", "nvm", "anms", "tnewt")
mychk <- opm(par=c(-1.2,1), fr, grr, method=mset)
cat("Summary output from opm\n")
print(summary(mychk))
cat("extract result for method 3\n")
print(opm2optimr(mychk, 3))
cat("Alternatively for method nvm\n")
print(opm2optimr(mychk, "nvm"))
cat("Bad inputs check for character method \n")
print(try(opm2optimr(mychk, "nvv")))
cat("Bad inputs check for out of range integer")
print(try(opm2optimr(mychk, 6)))
}
\keyword{nonlinear}
\keyword{optimize}
