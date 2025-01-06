\name{pd_check}
\alias{pd_check}
\encoding{UTF-8}
\title{Check Hessian matrix is positive definite by attempting a Cholesky decomposition.}
\description{We desire that Hessian matrices of objective functions for 
optimization are symmetric and positive definite. Symmetry is easy to check, but
we caution that evaluation or approximation of a Hessian often gives a slightly
asymmetric, and we commonly average the elements that are supposed to be equal,
thereby symmetrizing the proposed Hessian. 

The Cholesky decomposition will fail for a symmetric matrix which is NOT positive definite.
}
\usage{
   pd_check(A, pivot=TRUE, tol=1.e-07) 
}
\arguments{
  \item{A}{A symmetric matrix}
  \item{pivot}{Use pivoting in the Cholesky decomposition attempt}
  \item{tol}{A tolerance for zero}
}
\value{
  The function returns a logical TRUE if the decomposition can be carried out,
  FALSE otherwise.
}
\examples{
  n <- 10
  A <- pracma::moler(n)
  system.time(print(pd_check(A))) # returns FALSE
  system.time(print(eigen(A)$value))
  n <- 100
  A <- pracma::moler(n)
  system.time(print(pd_check(A))) # returns FALSE
  system.time(print(eigen(A)$value))
}
\keyword{positive-definite}
\keyword{Cholesky}

