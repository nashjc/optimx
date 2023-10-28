# opm2optimr function to extract
# par, value, counts, convergence
# do we want full optimr form (2 element counts, message and hessian), even
#  if we make those emtpy but they match structure?
#
# optim output -- minor modification
# par --    The best set of parameters found.
# value	--   The value of fn corresponding to par.
# counts	-- A two-element integer vector giving the number of calls to fn and gr 
#            respectively. This excludes those calls needed to compute the Hessian
#             even though the opm() result will have these counts
# convergence	-- An integer code. 0 indicates successful completion
# message	-- A character string which for optim() or optimr() may give additional 
#            information returned by the optimizer, or NULL.
#            Here will be "Result of conversion from opm() result"
# 
opm2optimr <- function(opmobj, rid){
# opmobj is the opm() output object (NOT the summary() result)
# rid is the row id. It can be name of solver in quotes, or an integer.
    nropm<-dim(opmobj)[1]
    if(is.numeric(rid)){
       if ((rid < 1) || (rid > nropm)) stop("invalid numeric rid=",rid)
       crid <- rownames(opmobj)[rid]
    }
    if(is.character(rid)){
       if (rid %in% rownames(opmobj)) {
          nrid <- which(rownames(opmobj) == rid)
       } else {
          stop("invalid character rid=",rid)
       }
    }
    row <- opmobj[rid,]
    # nparams + 8
    n <- length(row) - 8
    par <- row[1:n]
    value <- row[n+1]
    counts <- row[(n+2):(n+3)] # Just 2 elements
    ccode <- row[(n+5)]
    message <- paste("Result of conversion from opm() result for method ID=",rid,sep='')
    retval <- list(par=par, value=value, counts=counts, convergence=ccode, message=message)
    retval
}
