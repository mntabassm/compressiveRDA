#' crda.demo
#'
#' @title
#' Classification Example (demo) for Compressive Regularized Discriminant Analysis (CRDA) approach using the Ell2-RSCM estimator of the covariance matrix.
#' @aliases crda.demo
#'
#' @description
#' The function \code{crda.demo} performs classification using CRDA2 for one split of a real genomic dataset, Khan'2001.
#'
#' @return An object \code{res} of class \code{crda.demo} with the following attributes:
#' \item{funCall}{The call to the \code{crda.demo} function.}
#' \item{TER}{Test error rate (TER)}
#' \item{FSR}{Feature selection rate (FSR)}
#' \item{CT}{Average computational time (CT) in seconds.}

#' @author
#' Muhammad Naveed Tabassum and Esa Ollila, 2019.
#'
#' @seealso \code{\link{crda}}, \code{\link{crda.auxfuns}}
#'
#' @export
#' @examples
#' crda.demo()
#' crda.demo(prior = 'estimated')


crda.demo <- function(prior = 'uniform'){

    cat(sprintf("CRDA's demo (example): classification using a real genomic dataset, Khan'%d. \n", 2001))
    p <- dim(X)[1]
    n <- dim(X)[2]
    nt <- length(yt)

    ptm <- proc.time()
    obj <- crda(X, y, Xt = Xt, prior = prior, centerX = FALSE)
    CT <- as.numeric((proc.time()-ptm)[3])
    FSR <- obj$K / p
    TER <- sum(obj$yhat != yt) / nt	# Test Error Rate
    print(sum(obj$yhat))
    res <- list(funCall = match.call(), class = "crda.demo", TER = TER, FSR = FSR, CT = CT )
    return(res)
}
