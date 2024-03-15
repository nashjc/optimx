############### grnd.R ####################
grpracma <- function(par, userfn, ...) { # using grad from numDeriv
   tryg<-pracma::grad(userfn, par, ...)
}
############### end grnd ####################

