#' crda
#'
#' @title
#' Compressive Regularized Discriminant Analysis (CRDA) Method
#' @aliases crda
#'
#' @description
#' The function \code{crda} implements Compressive Regularized Discriminant Analysis (CRDA) approach and performs simultaneous feature selection and classification of high-dimensional data.
#' CRDA approach aims to address three facets of high-dimensional classification: namely, accuracy, computational complexity, and interpretability.
#'
#' @param X Training dataset, a pxn matrix with n-samples each having p-features.
#' @param y Labels for training dataset, an nx1 vector of whole numbers.
#' @param q Type of Lq,1 norm, default is Linf-norm.
#' @param al Regularization parameter.
#' @param K Joint-sparsity level.
#' @param Xt Test dataset, a pxnt matrix with nt-samples each of p-features.
#' @param prior Type of prior class probabilities, either 'uniform' (default) or 'estimated'.
#' @param centerX Flag for grand-mean centering of test dataset using grand-mean of training dataset.
#'
#' @return An object \code{obj} of class \code{crda} with the following attributes:
#' \item{funCall}{The call to the \code{crda} function.}
#' \item{prior}{Prior class probabilities.}
#' \item{varSelRate}{Feature selection rate (FSR).}
#' \item{selVarPos}{Position (i.e., index) of selected features.}
#' \item{coefMat}{Coefficient matrix before feature selection.}
#' \item{shrunkenCoefMat}{Shrunken (rowsparse) coefficient matrix.}
#' \item{const}{The constant part of discriminant function for CRDA method.}
#' \item{predTrainLabels}{Predicted labels for training dataset.}
#' \item{predTestLabels}{Optional: Predicted labels for test dataset, if it is available.}
#' \item{regparam}{Optional: The value of regularization parameter.}
#' \item{muX}{Optional: Grand-mean, i.e., row (feature) wise mean of training dataset.}
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2018
#'
#' @seealso \code{\link{crda.regparam}}, \code{\link{crda.cv}}
#'
#' @export
#' @examples
#' crda(X,y)
#' crda(X,y, Xt = testdata)
#' crda(X,y, q = 1, Xt = testdata)
#' crda(X,y, q = 2, Xt = Xt, centerX = TRUE)
#' crda(X,y, Xt = Xt, prior = 'estimated', centerX = TRUE)

crda <- function(X, y, q = Inf, al = NULL, K = NULL, Xt = NULL,
                 prior = 'uniform', centerX = FALSE){

  obj <- list(funCall = match.call(), class = "crda")

  if (is.null(al) || missing(al)) {
    al <- crda.regparam(X)    # Estimate of CRDA's regularization parameter using training dataset
    obj$regparam <- al
  }

  if (centerX) {
    muX <- apply(X, 1, mean)  # Grand-mean, i.e., row (feature) wise mean of training dataset, which will
    X <- X - muX              # be used later for grand-mean centering of test dataset when centerX=TRUE
    obj$muX <- muX
  }

  p <- dim(X)[1]
  n <- dim(X)[2]
  Y <- model.matrix( ~ factor(y)-1)
  M <- scale(X %*% Y, FALSE, table(y)) # Group-means matrix M

  if (prior == 'uniform') {
    prior <- rep(1/ncol(Y), ncol(Y))  # Uniform priors
  }else{
    prior <- table(y)/length(y)       # Estimated priors are relative sample proportions
  }
  obj$prior = prior

  Xc <- X - M[, unclass(factor(y))]    # Group-mean centering of training dataset
  xtx <- t(Xc) %*% Xc
  s <- svd(xtx)
  Dt2 <- s$d
  svalPos <- seq(Dt2)[Dt2 > 1e-5]      # Nonzero singular values' positions
  m <- length(svalPos)                 # Number of nonzero singular values OR rank of X
  D <- sqrt(Dt2[svalPos])
  U <- scale(Xc %*% s$u[, svalPos], FALSE, D) # As U = XVD^-1
  UtM <- t(U) %*% M

  eta <- (1-al) * ( sum(D^2)/(n*p) ) # (1-\alpha)*\eta
  D <- D / sqrt(n)
  regD <- 1/(al*D^2 + eta) -  1/eta
  regD <- as.vector(regD) * UtM

  B <- U %*% regD + M/eta            # Coefficient matrix B
  if (!is.finite(q)) {
    len <- apply(abs(B),1,max)          # Default: L_inf norm
  }else{
    len <- apply(abs(B)^q,1,sum)^(1/q)  # Lq norm
  }

  if(is.null(K) || missing(K)){
    K <- length(which(len > mean(len)))
  }
  I <- order(len, decreasing = TRUE)
  sflist <- I[1:K]
  Bsh <- B
  Bsh[setdiff(1:p,sflist), ] <- 0               # Setting (p-k)-features as zeros to get k-rowsparsity
  const <- -0.5*diag(t(M) %*% Bsh) + log(prior)	# 1xG vector of constant in RDA's Disc.Fn

  obj$varSelRate <- K/p
  obj$selVarPos <- sflist
  obj$coefMat <- B
  obj$shrunkenCoefMat <- Bsh
  obj$const <- const

  discFun <- scale(t(X) %*% Bsh, -const, FALSE)	# Predicted class posterior probabilities for training set
  obj$predTrainLabels <- max.col(discFun)

  if(!(is.null(Xt) || missing(Xt))){
    if (centerX) {
      Xt <- Xt - muX  # Grand-mean centering of test dataset using grand-mean of training dataset
    }
    discFun <- scale(t(Xt) %*% Bsh, -const, FALSE)  # Predicted class posterior probabilities for test set
    obj$predTestLabels <- max.col(discFun)
  }

  return(obj)
}
