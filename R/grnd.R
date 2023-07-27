############### grnd.R ####################
grnd<-function(par, userfn, ...) { # using grad from numDeriv
# with Richardson method
   tryg<-numDeriv::grad(userfn, par, ...)
}
############### end grnd ####################

