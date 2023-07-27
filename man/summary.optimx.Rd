\name{summary.optimx}
\alias{summary.optimx}
\encoding{UTF-8}
\title{Summarize optimx object}
\concept{minimization}
\concept{maximization}
\description{
  Summarize an \code{"optimx"} object.
}
\usage{
  \method{summary}{optimx}(object, order = NULL, par.select = TRUE, ...) 

}
\arguments{
 \item{object}{Object returned by \code{optimx}.}
 \item{order}{A column name, character vector of columns names, R 
    expression in terms of column names or a list of R expressions in terms of
    column names.   \code{NULL}, the default, means no re-ordering.
    \code{rownames} can be used to 
    alphabetic ordering by method name.
    \code{NULL}, the default, causes
    it not to be reordered.  Note that if 
    \code{follow.on} is TRUE re-ordering likely makes no sense. 
    The result is ordered by the \code{order}
    specification, each specified column in ascending order (except for 
    \code{value} which is in descending order if the optimization problem is a 
    maximization problem).}
 \item{par.select}{a numeric, character or logical vector selecting
   those \code{par} values to display. For example, \code{par=1:5} means
   display only the first 5 parameters.  Recycled so \code{par.select=FALSE}
   selects no parameters.}
 \item{\dots}{Further arguments to be passed to the function. Currently not
used.}
}
\details{
If \code{order} is specified then the result is reordered by the specified
columns, each in ascending order (except possibly for the
\code{value} column which is re-ordered in descending order for
maximization problems).}
\value{
  \code{summary.optimx} returns \code{object} with the rows ordered according
to \code{order} and with those parameters selected by \code{par.select}.
}
\examples{
ans <- optimx(fn = function(x) sum(x*x), par = 1:2)

# order by method name.
summary(ans, order = rownames)

# order by objective value. Do not show parameter values.
summary(ans, order = value, par.select = FALSE)

# order by objective value and then number of function evaluations
# such that objectives that are the same to 3 decimals are 
# considered the same.  Show only first parameter.
summary(ans, order = list(round(value, 3), fevals), par.select = 1)
}
\keyword{nonlinear}
\keyword{optimize}

