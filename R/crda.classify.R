#' crda.classify
#'
#' @title
#' Classification example by CRDA
#' @aliases crda.classify
#'
#' @description
#' The function \code{crda.cv} performs Classification by CRDA-variants
#'
#' @param L Number of training and test splits
#' @param nK Number of candidates in CV grid for joint-sparsity level
#' @param q Type of Lq,1 norm
#' @param prior Type of prior class probabilities
#' @param centerX Flag for centering the training dataset
#' @param comparison Flag for comparison with other methods
#'
#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2018
#'
#' @seealso \code{\link{crda}}, \code{\link{crda.setup3}}, \code{\link{crda.regparam}}, \code{\link{crda.cv}}
#'
#' @export
#' @examples
#' crda.classify()

crda.classify <- function(L = 10, nK = 10, q = Inf,
                          prior = 'uniform', centerX = FALSE, comparison = FALSE){

  if (!comparison) {
    algosList <- c("CRDA(Kub)","CRDA(Kcv)")
  }else{
    algosList <- c("CRDA(Kub)","CRDA(Kcv)","logit-ASPLS","SPCALDA","PLDA",
                   "SCRDA","NSC","DQDA","DLDA","varSelRF","lin-SVM")
  }
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
