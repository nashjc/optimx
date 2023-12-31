hesschk <- function(xpar, ffn, ggr, hhess, trace=0, testtol=(.Machine$double.eps)^(1/3), ...) {
   # check Hessian code in hess for function in fn and gradient in gr
   asymtol<-testtol # May want to change this.
   hessOK<-FALSE
   if (is.null(hhess)) {
      attr(hessOK, "nullhess")<-TRUE
      attr(hessOK, "asym")<-NA
      attr(hessOK, "ha")<-NA
      attr(hessOK, "hn")<-NA
      msg<-"Analytic Hessian not made available."
      attr(hessOK, "msg")<-msg
      if (trace > 0) cat(msg,"\n")
   } else {
      attr(hessOK, "nullhess")<-FALSE
      hname <- deparse(substitute(hhess))
      if (trace > 0) cat("Analytic hessian from function ",hname,"\n\n")
      ha <- hhess(xpar, ... ) # analytic hessian 
#      if (attr(ha,"inadmissible")) {
#         msg<-"Analytic Hessian inadmissible."
#         attr(hessOK, "msg")<-msg
#         if (trace > 0) cat(msg,"\n")
#      }
      if (is.null(ggr)) {
         ffn1 <- function(xpar) ffn(xpar, ...) # avoid dotarg name clash
         hn <- hessian(func=ffn1, x=xpar) 
      } else {
         ggr1<- function(xpar) ggr(xpar, ...) # avoid dotarg name clash
         hn <- jacobian(func=ggr1, x=xpar)
      }
      asym<-0.0 # to ensure defined
      if (!isSymmetric(hn)) {
          asym <- sum(abs(t(hn) - hn))/sum(abs(hn))
          asw <- paste("hn from hess() is reported non-symmetric with asymmetry ratio ", 
                  asym, sep = "")
          if (trace > 0) cat(asw, "\n")
          else warning(asw)
          if (asym > asymtol) {
             msg<-"Analytic Hessian not symmetric."
             attr(hessOK, "msg")<-msg
             if (trace > 0) cat(msg,"\n")
          } else hessOK <- TRUE
          hn <- 0.5 * (t(hn) + hn)
      }  # end if ! isSymmetric
      # Now test for equality. Again have to consider tolerance values
      fval<-ffn(xpar, ...) #  Could consider providing this externally
      if (max(abs(hn-ha))/(1 + abs(fval)) >= testtol) {
         hessOK<-FALSE
         msg<-paste("Analytic Hessian and numeric Hessian differ more than ", testtol,"")
         attr(hessOK, "msg")<-msg
         if (trace > 0) cat(msg,"\n")
      }
      attr(hessOK, "asym")<-asym
      attr(hessOK, "ha")<-ha
      attr(hessOK, "hn")<-hn
   } # end gradient/hessian tests
   hessOK
}

