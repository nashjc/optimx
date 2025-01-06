\name{optimr2opm}
\alias{optimr2opm}
\encoding{UTF-8}
\title{Add a single optimr() solution to a opm() result set}
\concept{minimization}
\concept{maximization}
\description{
  A function that attempts to add a row to an opm() result matrix.
  opm() is a multi-method optimization wrapper.
}
\usage{
optimr2opm(ans, opmmat)
}
\arguments{
 \item{ans}{A result object from the optimr() function for a particular
 \code{method}.}
 \item{opmmat}{the object returned by opm() from trying to solve an
  optimization problem with multiple solvers}
}
\value{
  An object in the form of an opm() multiple solver result with either one
  row if \code{opmobj} does not exist, or one new row at the bottom if it
  does.
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
mset <- c("ncg", "nvm", "anms")
x0<-c(-1.2,1)
mychk <- opm(par=x0, fr, grr, method=mset)
cat("Summary output from opm\n")
print(summary(mychk))
myans<-optimr(x0, fr, grr, method="tnewt")
cat("Add the solution from method 'tnewt'\n")
mychk2<-optimr2opm(myans, mychk)
print(summary(mychk2))
}
\keyword{nonlinear}
\keyword{optimize}
