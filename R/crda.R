#' crda
#'
#' @title
#' Compressive regularized discriminant analysis (CRDA) approach based CRDA2 method
#' @aliases crda
#'
#' @description
#' CRDA approach performs simultaneous feature selection and classification of high-dimensional data.
#' The function \code{crda} classifies each column of the test data set Xt (p x N) into one of
#' the G classes. Test data set Xt and training data set X must have the
#' same number of rows (features or variables). Vector y is a
#' class variable of training data. Its unique values define classes; each
#' element defines the class to which the corresponding column of X belongs.
#' The input y is a numeric vector with integer elements ranging from
#' 1,2,..,G, where G is the number of classes. Note that y must have the
#' same number of rows as there are columns in X. The output yhat indicates
#' the class to which each column of Xt has been assigned. Also yhat is Nx1
#' vector of integers ranging from 1,2,...,G.
#' CRDA approach aims to address three facets of high-dimensional classification: namely, accuracy, computational complexity, and interpretability.
#'
#' @param X Training dataset, a pxn matrix with n-samples each having p-features.
#' @param y Labels for training dataset, an nx1 vector of whole numbers.
#' @param Xt Test dataset, a pxnt matrix with nt-samples each of p-features.
#' @param q scalar (>=1) or string 'var' or 'inf' or 'cv'. If q is a real scalar, then it must be >= 1 and it
#' denotes the L_q-norm to be used in the hard thresholding operator H_K(B,phi) in the CRDA method.
#' If q is equal to a  string 'inf' (or 'var') then L_infty norm will be used (or sample variance) will be  used
#' as the hard-thresholding  selector function. If q is equal to a string 'cv', then CV is used to select the best
#' hard-thresholding selector function  among the L_1-, L_2-, L_inf-norm and the sample variance.
#' @param prior Type of prior class probabilities, either 'uniform' (default) or 'estimated'.
#' @param kgrid A grid of candidate values for joint-sparsity level.
#' @param nK Number of candidate values in the grid of joint-sparsity level.
#' @param nfolds Number of folds for cross-validation scheme.
#' @param centerX Flag for grand-mean centering of test dataset using grand-mean of training dataset.
#'
#' @return An object \code{obj} of class \code{crda} with the following attributes:
#' \item{funCall}{The call to the \code{crda} function.}
#' \item{prior}{Prior class probabilities.}
#' \item{B}{Coefficient matrix before feature selection.}
#' \item{M}{Group-means matrix.}
#' \item{Ind}{Indices of the length of the row vectors of B organized in descending order, where the length is determined by optional argument 'q'.}
#' \item{K}{Joint-sparsity level.}
#' \item{yhat}{Predicted labels for test dataset}
#' \item{muX}{Optional: Grand-mean, i.e., row (feature) wise mean of training dataset.}
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2019.
#'
#' @seealso \code{\link{crda.auxfuns}}, \code{\link{crda.demo}}
#'
#' @importFrom caret createFolds
#' @export
#' @examples
#' crda(X,y)
#' crda(X,y, Xt = testdata)
#' crda(X,y, Xt = testdata, q = 1)
#' crda(X,y, Xt = Xt, q = 2, centerX = FALSE)
#' crda(X,y, Xt = Xt, prior = 'estimated', centerX = TRUE)

crda <- function(X, y, Xt = NULL, q = 'cv', prior = 'uniform', kgrid = NULL,
                 nK = 10, nfolds = 5, centerX = TRUE){

  obj <- list(funCall = match.call(), class = "crda")
  cat(sprintf("CRDA using the Ell%d-RSCM estimator of the covariance matrix. \n", 2))

  if (centerX) {
    muX <- apply(X, 1, mean)  # Grand-mean, i.e., row (feature) wise mean of training dataset, which will
    X <- X - muX              # be used later for grand-mean centering of test dataset when centerX=TRUE
    obj$muX <- muX
  }

  p <- dim(X)[1]
  n <- dim(X)[2]

  G <- base::length(unique(y))
  if (prior == 'uniform') {
    prior <- rep(1/G, G)  # Uniform priors
  }else{
    prior <- table(y)/base::length(y)       # Estimated priors are relative sample proportions
  }
  obj$prior = prior

  obj0 <- CRDA_coefmat(X, y)
  B <- obj0$B
  obj$B <- B
  M <- obj0$M
  obj$M <- M

  if (q>=1 && is.numeric(q)) {
    len <- cbind(apply(abs(B)^q,1,sum)^(1/q))
  }else{
    len <- cbind(switch(q, "var" = apply(B,1,var), "inf" = apply(abs(B),1,max), "cv" = NULL))
  }
  if (is.null(len)){
    len <- cbind(apply(B,1,var),  apply(abs(B),1,max), apply(abs(B),1,sum), sqrt(apply(abs(B)^2,1,sum)))
    qvals <- list("inf", "var", 2,1)
  }else{
    qvals <- list(q)
  }
  Kup <- apply(len>apply(len,2,mean),2,sum)
  Kup <- min(Kup)
  if (is.null(kgrid)){
    kgrid <- unique(round(exp(seq(log(0.05*p),log(Kup), length.out = nK))))
  }

  flds  <- caret::createFolds(y, k = nfolds, list = TRUE, returnTrain = FALSE)
  cverr <- matrix(0, nK, dim(len)[2])
  Ind <- matrix(0, p, dim(len)[2])

  for (fnum in 1:nfolds) {

    testIndx <- flds[[fnum]]
    trainIndx <- setdiff(1:base::length(y),testIndx)
    Xii <- X[,trainIndx]
    yii <- y[trainIndx]
    Xho <- cbind(X[,testIndx])
    yhat <- matrix(0, dim(Xho)[2], dim(len)[2])

    for (kk in 1:nK) {

      K = kgrid[kk]

      if (kk!=1) {

        if (q!='cv') {
          res <- CRDA0(Xho,Xii,yii,K,q,prior,Bii,Mii,Ind)
          yhat <- res$yhat
          }else{
            res <- CRDA0(Xho,Xii,yii,K,qvals[[1]],prior,Bii,Mii,Ind[,1])
            yhat[,1] <- res$yhat
            res <- CRDA0(Xho,Xii,yii,K,qvals[[2]],prior,Bii,Mii,Ind[,2])
            yhat[,2] <- res$yhat
            res <- CRDA0(Xho,Xii,yii,K,qvals[[3]],prior,Bii,Mii,Ind[,3])
            yhat[,3] <- res$yhat
            res <- CRDA0(Xho,Xii,yii,K,qvals[[4]],prior,Bii,Mii,Ind[,4])
            yhat[,4] <- res$yhat
          }

      }else{
        if (q!='cv') {
          res <- CRDA0(Xho,Xii,yii,K,q,prior)
          yhat <- res$yhat
          Bii <- res$B
          Mii <- res$M
          Ind <- res$Ind
        }else{
          res <- CRDA0(Xho,Xii,yii,K,qvals[[1]],prior)
          Bii <- res$B
          Mii <- res$M
          yhat[,1] <- res$yhat
          Ind[,1]  <- res$Ind
          res <- CRDA0(Xho,Xii,yii,K,qvals[[2]],prior,Bii,Mii)
          yhat[,2] <- res$yhat
          Ind[,2]  <- res$Ind
          res <- CRDA0(Xho,Xii,yii,K,qvals[[3]],prior,Bii,Mii)
          yhat[,3] <- res$yhat
          Ind[,3]  <- res$Ind
          res <- CRDA0(Xho,Xii,yii,K,qvals[[4]],prior,Bii,Mii)
          yhat[,4] <- res$yhat
          Ind[,4]  <- res$Ind
        }
        }

      cverr[kk,] <- cverr[kk,] + apply(yhat != y[testIndx],2,sum)  / base::length(testIndx)
    }
  }

  cverr = cverr/nfolds
  idxmin <- which(cverr <= min(cverr), arr.ind = TRUE)
  idxK <- idxmin[,1]
  idxq <- idxmin[,2]
  idxK0 = min(idxK)
  # if there are many q values that had small CV err, then
  # determine the best as the one having small mean CV err
  if (sum(idxK==idxK0) != 1) {
    bst_indx = idxq[idxK==idxK0]
    cverr0 <- apply(cverr[,bst_indx],2,mean)
    idxq0 <- bst_indx[ which(cverr0 <= min(cverr0), arr.ind = TRUE)[1] ]
  } else{
        idxq0 = idxq[idxK==idxK0]
  }
  K <- kgrid[idxK0]
  res <- CRDA0(Xt,X,y,K,qvals[[idxq0]],prior,B,M)
  obj$Ind  <- res$Ind
  obj$K  <- res$K
  obj$yhat <- res$yhat

  return(obj)
}

