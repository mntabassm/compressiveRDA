#' crda.classify
#'
#' @title
#' Classification Example using Compressive Regularized Discriminant Analysis (CRDA) Method
#' @aliases crda.classify
#'
#' @description
#' The function \code{crda.cv} performs classification using CRDA-variants for a partially synthetic dataset.
#'
#' @param L Number of runs, i.e., training and test splits.
#' @param nK Number of candidates in 5-fold CV-grid for finding joint-sparsity level.
#' @param q Type of Lq,1 norm, default is Linf-norm.
#' @param prior Type of prior class probabilities, either 'uniform' (default) or 'estimated'.
#' @param centerX Flag for grand-mean centering of test dataset using grand-mean of training dataset.
#'
#' @return An object \code{res} of class \code{crda.classify} with the following attributes:
#' \item{funCall}{The call to the \code{crda.classify} function.}
#' \item{ACT}{Average computational time (ACT) in seconds over L-runs.}
#' \item{AveTER}{Average test error rate (TER) over L-runs.}
#' \item{AveFSR}{Average feature selection rate (FSR) over L-runs.}
#' \item{AveFPR}{Average false positive rate (FPR) over L-runs.}
#' \item{AveFNR}{Average false negative rate (FNR) over L-runs.}

#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2018
#'
#' @seealso \code{\link{crda}}, \code{\link{crda.setup3}}, \code{\link{crda.regparam}}, \code{\link{crda.cv}}
#'
#' @export
#' @examples
#' crda.classify()
#' crda.classify(L = 1)
#' crda.classify(L = 1, nK = 5)
#' crda.classify(L = 1, prior = 'estimated')
#' crda.classify(L = 1, centerX = TRUE)


crda.classify <- function(L = 10, nK = 10, q = Inf,
                          prior = 'uniform', centerX = FALSE){

  algosList <- c("CRDA(Kub)","CRDA(Kcv)")
  CT <- FSR <- TER <- FPR <- FNR <- matrix(NaN, L, length(algosList))
  TERnaive <- rep(NaN, L)

  set.seed(1)
  for (mc in 1:L){

    cat(sprintf("Run= %d | Generating training and test datasets... \n", mc))
    dataset <- crda.setup3()
    X  <- dataset$X
    y  <- dataset$y
    yt <- dataset$yt
    p <- dim(X)[1]
    n <- dim(X)[2]
    nt <- length(yt)
    Tp <- length(dataset$truePos)
    G <- length(unique(y))
    flds <- caret::createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

    opt_al <- crda.regparam(X) # estimating regularization parameter using training dataset

    cat(sprintf("Classification by CRDA(Kub)... "))
    algo <- 1
    ptm <- proc.time()
    out <- crda(X, y, q = q, al = opt_al, Xt = dataset$Xt, prior = prior, centerX = centerX)
    CT[mc,algo] <- as.numeric((proc.time()-ptm)[3])
    FSR[mc,algo] <- out$varSelRate
    TER[mc,algo] <- sum(out$predTestLabels != yt) / nt	# Test Error Rate
    Kub <- length(out$selVarPos)
    C <- sum(dataset$truePos %in% out$selVarPos)
    FPR[mc,algo] <- (Tp-C) / Tp
    FNR[mc,algo] <- (Kub-C) / (p-Tp)
    cat(sprintf("Done! | "))

    if(Kub>0.1*p){
      Kgrid <- round(exp(seq(log(0.05*p),log(Kub), length.out = nK)))
    }else{
      Kgrid <- round(exp(seq(log(0.05*p),log(0.4*p), length.out = nK)))
    }

    cat(sprintf("Classification by CRDA(Kcv)... "))
    algo <- 2
    ptm <- proc.time()
    outCV <- crda.cv(X, y, q = q, al = opt_al, flds = flds, Kgrid = Kgrid, prior = prior, centerX = centerX)
    out <- crda(X, y, q = q, al = opt_al, K = outCV$K, Xt = dataset$Xt, prior = prior, centerX = centerX)
    CT[mc,algo] <- as.numeric((proc.time()-ptm)[3])
    FSR[mc,algo] <- out$varSelRate
    TER[mc,algo] <- sum(out$predTestLabels != yt) / nt	# Test Error Rate
    Kcv <- length(out$selVarPos)
    C <- sum(dataset$truePos %in% out$selVarPos)
    FPR[mc,algo] <- (Tp-C) / Tp
    FNR[mc,algo] <- (Kcv-C) / (p-Tp)
    cat(sprintf("Done! \n"))

    ## Baseline 'naive' method, which classifies all of the test observations to the most frequent class in the training set.
    TERnaive[mc] <- sum(rep(which.max(as.numeric(table(y))), length(yt)) != yt) / nt
  }

  res <- list(funCall = match.call(), class = "crda.classify", ACT = apply(CT, 2, mean),
              AveTER = apply(TER, 2, mean),  AveFSR = apply(FSR, 2, mean),
              AveFPR = apply(FPR, 2, mean),  AveFNR = apply(FNR, 2, mean) )
  return(res)
}
