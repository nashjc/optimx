\name{scalechk}
\alias{scalechk}
\encoding{UTF-8}
\title{Check the scale of the initial parameters and bounds input to an optimization code
   used in nonlinear optimization}
   
\description{
Nonlinear optimization problems often have different scale for different
parameters. This function is intended to explore the differences in scale. It is, however,
an imperfect and heuristic tool, and could be improved.

At this time scalechk ignores parameters and bounds for fixed (masked) parameters
for calculations of scaling measures. The rationale for this is that such parameters
are outside the optimization process.
}

\usage{
   scalechk(par, lower = lower, upper = upper, bdmsk=NULL, dowarn = TRUE)
}

\arguments{
 \item{par}{A numeric vector of starting values of the optimization function parameters.}
 \item{lower}{A vector of lower bounds on the parameters.}
 \item{upper}{A vector of upper bounds on the parameters.}
 \item{bdmsk}{An indicator vector, having 1 for each parameter that is "free" or
     unconstrained, and 0 for any parameter that is fixed or MASKED for the
     duration of the optimization. May be NULL.}
 \item{dowarn}{Set TRUE to issue warnings. Othwerwise this is a silent routine.
               Default TRUE.}
}
\details{
     The scalechk function will check that the bounds exist and are admissible, 
     that is, that there are no lower bounds that exceed upper bounds. 

     NOTE: Free paramters outside bounds are adjusted to the nearest bound.
     We then set parchanged = TRUE which implies the original parameters
     were infeasible.

}

\value{
  A list with components:

   # Returns:
   #   list(lpratio, lbratio) -- the log of the ratio of largest to smallest parameters
   #      and bounds intervals (upper-lower) in absolute value (ignoring Inf, NULL, NA)

  \item{lpratio}{The log of the ratio of largest to smallest parameters
       in absolute value (ignoring Inf, NULL, NA)}
  \item{lbratio}{The log of the ratio of largest to smallest bounds intervals 
        (upper-lower) in absolute value (ignoring Inf, NULL, NA)}
}
\examples{
#####################
  par <- c(-1.2, 1)
  lower <- c(-2, 0)
  upper <- c(100000, 10)
  srat<-scalechk(par, lower, upper,dowarn=TRUE)
  print(srat)
  sratv<-c(srat$lpratio, srat$lbratio)
  if (max(sratv,na.rm=TRUE) > 3) { # scaletol from ctrldefault in optimx
     warnstr<-"Parameters or bounds appear to have different scalings.\n
     This can cause poor performance in optimization. \n
     It is important for derivative free methods like BOBYQA, UOBYQA, NEWUOA."
     cat(warnstr,"\n")
  }

}
\keyword{nonlinear}
\keyword{optimize}
\keyword{upper}
\keyword{lower}
\keyword{bound}
\keyword{mask}
