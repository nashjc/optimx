\name{checksolver}
\alias{checksolver}
\alias{checkallsolvers}
\encoding{UTF-8}
\title{Test if requested solver is present}
\concept{minimization}
\concept{maximization}
\description{
  Test if requested solver is present.
}
\usage{
   checksolver(method, allmeth, allpkg)
   checkallsolvers()
}
\arguments{
 \item{method}{Character string giving the name of the solver requested.}
 \item{allmeth}{Character vector giving the names of the methods optimr can use.}
 \item{allpkg}{Character vector giving the names of the packages where
     solvers are found.}
}
\value{
  If the solver defined by character string in \code{method} is
   available, then \code{checksolver} returns this string, else NULL.
   
  \code{checkallsolvers()} returns a vector of strings that are the
  names of missing solvers, else NULL if all solvers specified in \code{allmeth}
  are present where \code{allmeth} is returned from a call to \code{ctrldefault(n)} 
  where \code{n} is some integer.
}
\examples{
   allmeth <- c("Rvmmin", "nlminb","ipopttest")
   allpkg <- c("Rvmmin", "stats","ipoptr")
   
   print(checksolver("nlminb", allmeth, allpkg))
   # If Rvmmin NOT available, get msg that PACKAGE not available.
   print(checksolver("Rvmmin", allmeth, allpkg))
   # Get message that SOLVER not found
   print(checksolver("notasolver", allmeth, allpkg))

}
\keyword{nonlinear}
\keyword{optimize}

