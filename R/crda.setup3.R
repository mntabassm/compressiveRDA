
#' crda.setup3
#'
#' @aliases crda.setup3
#'
#' @description
#' This function generates partially synthetic data for setup 3 of the CRDA-paper.
#'
#' @return X,y, Xt,yt, truepos.
#' @importFrom mvtnorm rmvnorm
#' @examples
#' crda.setup3()

crda.setup3 <- function() {

  # load('khan2001.RData')
  # data(khan2001)
  ndf <- 100
  truePos <- 1:ndf
  p <- ncol(Xo)-ndf
  Xo[ ,(1:p)+ndf] <- mvtnorm::rmvnorm(n = nrow(Xo), rep(0,p), diag(1e-2,p))
  n <- length(yo)

  smp_siz = ceiling(0.6*n)  # creates a value for dividing the data into train and test. In this case the value is defined as 75% of the number of rows in the dataset
  all_labels <- FALSE
  while (!all_labels) {
    train_ind = sample(seq_len(n),size = smp_siz)  # Randomly identifies therows equal to sample size ( defined in previous instruction) from  all the rows of Smarket dataset and stores the row number in train_ind
    y <- yo[train_ind]
    yt <- yo[-train_ind]
    all_labels = setequal( unique(yo), Reduce(intersect, list(y,yt)) )
  }
  X  = Xo[train_ind, ]    # creates the training dataset with row numbers stored in train_ind
  Xt = Xo[-train_ind, ]  # creates the test dataset excluding the row numbers mentioned in train_ind
  dataset <- list(X = t(X), y = y, Xt = t(Xt), yt = yt, truePos = truePos)
  return(dataset)
}
