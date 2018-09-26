#' crda.cv
#'
#' @title
#' Cross-validation for CRDA
#' @aliases crda.cv
#'
#' @description
#' The function \code{crda.cv} performs cross-validation
#'
#' @param X Training dataset
#' @param y Labels for training dataset
#' @param q Type of Lq,1 norm
#' @param al regularization parameter
#' @param prior Type of prior class probabilities
#' @param flds Folds for cross-validation
#' @param nFolds number of folds
#' @param Kgrid A grid having candidate values of joint-sparsity parameter
#' @param centerX Flag for centering the training dataset
#' @param plotCV Flag for plotting the results
#'
#' @return An object of class \code{crda.cv} with the following attributes
#' \item{funCall}{The call for CRDA}
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2018
#'
#' @seealso \code{\link{crda.regparam}}, \code{\link{crda}}
#'
#' @importFrom caret createFolds
#' @examples
#' crda.cv(X,y) when using defaults
#' crda.cv(X,y, q=1,nFolds=10,centerX=TRUE)

crda.cv <- function(X, y, q = Inf, al = NULL, prior='uniform', flds = NULL, nFolds = 5,
                    Kgrid = round(exp(seq(log(0.05*nrow(X)),log(0.5*nrow(X)), length.out = 10))),
                    centerX = FALSE, plotCV = FALSE){

  if (is.null(al) || missing(al)) {
    al <- crda.regparam(X) # Estimate of CRDA's regularization parameter using training dataset
  }

  if (is.null(flds) || missing(flds)) {
    flds  <- caret::createFolds(y, k = nFolds, list = TRUE, returnTrain = FALSE)
  }else{
    nFolds <- length(flds)
  }
  nK    <- length(Kgrid)
  cvMER <- resubMER <- rep(0, nK)
  for (kk in 1:nK) {

    for (fnum in 1:nFolds) {
      testIndx <- flds[[fnum]]
      trainIndx <- setdiff(1:length(y),testIndx)
      out <- crda(X[,trainIndx],y[trainIndx], q = q, al = al, K = Kgrid[kk],
                  Xt=X[,testIndx], prior=prior, centerX=centerX)

      cvMER[kk]    <- cvMER[kk]    + sum(out$predTestLabels  != y[testIndx])  / length(testIndx)
      resubMER[kk] <- resubMER[kk] + sum(out$predTrainLabels != y[trainIndx]) / length(trainIndx)
    }
  }
  Kcv <- Kgrid[which.min(cvMER)]
  obj <- list(funCall = match.call(), Kgrid = Kgrid, folds = flds,
              cvMER = cvMER, resubMER = resubMER, regparam = al, K = Kcv, class = "crda.cv")
  return(obj)

  if(plotCV){
    lwd <- 1.5
    plotchar <- c(NA,NA,15)
    colors <- c("blue","black")
    FSR <- Kgrid/nrow(X)
    xrange <- c(min(FSR), max(FSR))
    yrange <- c(0, max(max(cvMER), max(resubMER)))
    plot(xrange, yrange, type="n", xlab="Feature selection rate", ylab="Error rate" )
    lines(FSR, resubMER, type="l", lwd=lwd, lty=1, col=colors[1], pch=plotchar[1])
    lines(FSR, cvMER,    type="l", lwd=lwd, lty=2, col=colors[2], pch=plotchar[2])
    lines(FSR[which.min(cvMER)], min(cvMER), type="p", col=colors[1], pch=plotchar[3])
    legend("top", horiz = TRUE, inset = c(0,-0.1), xpd = TRUE, bty = 'n',
           legend=c("reSub-MER","CV-MER","optFSR"), cex=1, col=colors, lty=c(1,2,0),pch=plotchar)
  }
}
