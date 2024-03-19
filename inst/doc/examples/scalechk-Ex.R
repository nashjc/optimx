#####################
library(optimx)
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

# test infinite values
par<-c(Inf, -Inf, 1, 1000)
lower<- -Inf
upper<- Inf
srat2<-try(scalechk(par, lower, upper, dowarn=TRUE))
print(srat2)

