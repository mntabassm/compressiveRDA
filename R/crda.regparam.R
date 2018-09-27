#' crda.regparam
#'
#' @title
#' Estimation of Regularization Parameter.
#' @aliases crda.regparam
#'
#' @description
#' The function \code{crda.regparam} estimates the regularization parameter by a closed-form solution.
#'
#' @param X Training dataset, a pxn matrix with n-samples each having p-features.
#'
#' @return The value of regularization parameter.
#' \item{al}{regularization parameter}
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2018
#'
#' @seealso \code{\link{crda}}, \code{\link{crda.cv}}
#'
#' @export
#' @examples
#' crda.regparam(X)


crda.regparam <- function(X){

  cat(sprintf("Estimating reg param... "))
  p <- dim(X)[1]
  n <- dim(X)[2]
  muX <- apply(X, 1, mean)
  Xc <- X - muX

  xtx <- t(Xc) %*% Xc
  s <- svd(xtx)
  Dt2 <- s$d
  svalPos <- seq(Dt2)[Dt2 > 1e-6]
  m <- length(svalPos)
  D <- sqrt(Dt2[svalPos])
  tr = sum(D^4)/(sum(D^2)^2)

  kurt1n <- (n-1)/((n-2)*(n-3))
  vari <- apply(Xc^2, 1, mean)
  g2 <- apply(Xc^4, 1, mean) / (vari^2)-3
  g2[which(is.nan(g2))] <- 0              # since g2 can have NaN due to 0/0 devision
  G2 <- kurt1n*((n+1)*g2 + 6)
  kurtest <- mean(G2)
  kappahat <- (1/3)*kurtest
  ka_lb = -2/(p+2)
  ka_lb = ka_lb + (abs(ka_lb))/40
  kappahat <- max(ka_lb, kappahat)

  tau2 <- kappahat/n
  tau1 <- 1/(n-1) + tau2
  a = tau1 / (1+tau2)
  b = (1+tau2) / (1+tau1*(1-2*tau1)+tau2*(2+tau1+tau2))
  gam = (b*p)*(tr - a)

  gam1 <- min(p-1, max(0,gam-1))
  al = gam1 / (gam1 + (gam+p)*tau1 + gam*tau2)
  alEll2 = min(1, max(0,al))
  cat(sprintf("Done! | alpha(Ell2) is %f \n", alEll2))

  return(alEll2)
}
