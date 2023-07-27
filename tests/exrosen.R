# Extended Rosenbrock Function ex_rosen.R from 
# https://github.com/jlmelville/funconstrain by jlmelville
# Test function 21 from the More', Garbow and Hillstrom paper.
#
# The objective function is the sum of \code{m} functions, each of \code{n}
# parameters.
#
xrosn.f = function(par) {
      n <- length(par)
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }

      fsum <- 0
      for (i in 1:(n / 2)) {
        p2 <- 2 * i
        p1 <- p2 - 1

        f_p1 <- 10 * (par[p2] - par[p1] ^ 2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
      }

      fsum
}

xrosn.g = function(par) {
      n <- length(par)
      if (n %% 2 != 0) {
        stop("Extended Rosenbrock: n must be even")
      }

      grad <- rep(0, n)
      for (i in 1:(n / 2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]

        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]

        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
      }

      grad
}

