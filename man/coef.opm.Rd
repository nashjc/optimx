\name{coef}
\alias{coef<-}
\alias{coef.opm}
\alias{coef<-.opm}
\alias{coef.optimx}
\alias{coef<-.optimx}
\encoding{UTF-8}
\title{Summarize opm object}
\concept{minimization}
\concept{maximization}
\description{
  Summarize an \code{"opm"} object.
}
\usage{
  \method{coef}{opm}(object, ...) 
  \method{coef}{opm}(x) <- value 
}
\arguments{
 \item{object}{Object returned by \code{opm}.}
 \item{\dots}{Further arguments to be passed to the function. Currently not
   used.}
 \item{x}{An \code{opm} object.}
 \item{value}{Set parameters equal to this value.}
}
\value{
  \code{coef.opm} returns the best parameters found by each method that
  returned such parameters. The returned coefficients are in the form of a
  matrix with the rows named by the relevant methods and the columns named
  according to parameter names provided by the user in the vector of starting
  values, or else by "p1", "p2", ..., if names are not provided.
}
\examples{
ans <- opm(fn = function(x) sum(x*x), par = 1:2, method="ALL", control=list(trace=0))
print(coef(ans))

ansx <- optimx(fn = function(x) sum(x*x), par = 1:2, control=list(all.methods=TRUE, trace=0))
print(coef(ansx))


\dontrun{
proj <- function(x) x/sum(x)
f <- function(x) -prod(proj(x))
ans <- opm(1:2, f)
print(ans)
coef(ans) <- apply(coef(ans), 1, proj)
print(ans)
}

}
\keyword{nonlinear}
\keyword{optimize}

