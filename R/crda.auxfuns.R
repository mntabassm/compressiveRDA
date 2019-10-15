
#' crda.auxfuns
#'
#' @title
#' Auxiliary functions used by the CRDA approach.
#' @aliases crda.auxfuns
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2019.
#'
#' @seealso \code{\link{crda}}, \code{\link{crda.demo}}

# ---------------------------------------------------------
# CRDA_coefmat Computes the CRDA coefficinent matrix B = Sigma^-1 * M
CRDA_coefmat <- function(X, y, al = NULL, M = NULL){

  obj <- list(funCall = match.call(), class = "CRDA_coefmat")

  if (is.null(al) || missing(al)) {
    al <- crda.regparam(X)    # Estimate of CRDA's regularization parameter using training dataset
    obj$regparam <- al
  }

  p <- dim(X)[1]
  n <- dim(X)[2]
  if (is.null(M) || missing(M)){
    Y <- model.matrix( ~ factor(y)-1)
    M <- scale(X %*% Y, FALSE, table(y)) # Group-means matrix M
  }
  obj$M <- M

  Xc <- X - M[, unclass(factor(y))]    # Group-mean centering of training dataset
  xtx <- t(Xc) %*% Xc
  s <- svd(xtx)
  Dt2 <- s$d
  svalPos <- seq(Dt2)[Dt2 > 1e-5]      # Nonzero singular values' positions
  m <- base::length(svalPos)                 # Number of nonzero singular values OR rank of X
  D <- sqrt(Dt2[svalPos])
  U <- scale(Xc %*% s$u[, svalPos], FALSE, D) # As U = XVD^-1
  UtM <- t(U) %*% M
  eta <- (1-al) * ( sum(D^2)/(n*p) ) # (1-\alpha)*\eta
  D <- D / sqrt(n)
  regD <- 1/(al*D^2 + eta) -  1/eta
  regD <- as.vector(regD) * UtM
  obj$B <- U %*% regD + M/eta            # Coefficient matrix B

  return(obj)
}

# ---------------------------------------------------------
# crda.regparam estimates the regularization parameter by a closed-form solution, in the Ell2-RSCM estimator of the covariance matrix.
crda.regparam <- function(X, centerX = TRUE){

  # cat(sprintf("Estimating reg param... "))
  p <- dim(X)[1]
  n <- dim(X)[2]
  if (centerX) {
    muX <- apply(X, 1, mean)  # Grand-mean, i.e., row (feature) wise mean of training dataset, which will
    X <- X - muX              # be used later for grand-mean centering of test dataset when centerX=TRUE
  }

  xtx <- t(X) %*% X
  s <- svd(xtx)
  Dt2 <- s$d
  svalPos <- seq(Dt2)[Dt2 > 1e-6]
  m <- length(svalPos)
  D <- sqrt(Dt2[svalPos])
  tr = sum(D^4)/(sum(D^2)^2)

  kurt1n <- (n-1)/((n-2)*(n-3))
  vari <- apply(X^2, 1, mean)
  g2 <- apply(X^4, 1, mean) / (vari^2)-3
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
  # cat(sprintf("Done! | alpha(Ell2) is %f \n", alEll2))

  return(alEll2)
}

# ---------------------------------------------------------
# CRDA0 is the basic utility function needed by the CRDA function.
CRDA0 <- function(Xt = NULL, X, y, K = NULL, q = 'inf', prior = table(y)/base::length(y), B = NULL, M = NULL, Ind = NULL){

  obj <- list(funCall = match.call(), class = "CRDA0")

  if (is.null(B)){
    obj0 <- CRDA_coefmat(X, y)
    B <- obj0$B
    M <- obj0$M
  }
  obj1 <- hard_threshold(B,K,q,Ind)
  Best <- obj1$B

  const <- -0.5*diag(t(M) %*% Best) + log(prior)	# 1xG vector of constant in RDA's Disc.Fn
  dlda <- scale(t(Xt) %*% Best, -const, FALSE)	# Predicted class posterior probabilities for test set
  obj$yhat <- max.col(dlda)
  obj$B <- B
  obj$M <- M
  obj$K <- obj1$K
  obj$Ind <- obj1$Ind
  return(obj)
}

# ---------------------------------------------------------
# Hard-thresholding operator H_K(B,q) used in the compressive regularized (linear) discriminant analysis (CRDA).
hard_threshold <- function(B, K = NULL, q = 'inf', Ind = NULL){

  obj <- list(funCall = match.call(), class = "hard_threshold")

  if (q>=1 && is.numeric(q)) {
    len <- cbind(apply(abs(B)^q,1,sum)^(1/q))
  }else{
    len <- cbind(switch(q, "var" = apply(B,1,var), "inf" = apply(abs(B),1,max)))
  }
  if(is.null(K) || missing(K)){
    K <- base::length(which(len > mean(len)))
  }
  Ind <- order(len, decreasing = TRUE)
  obj$Ind <- Ind
  B[setdiff(1:dim(len)[1],Ind[1:K]), ] <- 0       # Setting (p-k)-features as zeros to get k-rowsparsity
  obj$B <- B
  obj$K <- K
  return(obj)
}


