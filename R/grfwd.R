grfwd <- function(par, userfn, fbase=NULL, env=optsp, ...) {
   # Forward different gradient approximation
   # ?20220223 -- helicalfn gives poor results as dx ~= par. Changed deps in zzz.R
   eps<-env$deps
   if (is.null(fbase)) fbase <- userfn(par, ...)  #  function value at par
   df <- rep(NA, length(par))
   teps <- eps * (abs(par) + eps)
   dx <- par
   for (i in 1:length(par)) {
      dx[i] <- dx[i] + teps[i]
      df[i] <- (userfn(dx, ...) - fbase)/teps[i]
   }
   df
}

