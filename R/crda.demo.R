#' crda.demo
#'
#' @title
#' A demo (example) on classification of a real genomic dataset
#' @aliases crda.demo
#'
#' @description
#' The function \code{crda.demo} performs classification using compressive regularized
#' discriminant analysis (CRDA) approach for one split of a real genomic dataset, Khan'2001.
#'
#' @return An object \code{res} of class \code{crda.demo} with the following attributes:
#' \item{funCall}{The call to the \code{crda.demo} function.}
#' \item{TER}{Test error rate (TER) in terms of percentage of misclassifications.}
#' \item{FSR}{Feature selection rate (FSR) in terms of percentage of features used.}
#' \item{CT}{Computational time (CT) in seconds.}
#'
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
    FSR <- round(100*obj$K / p, digits = 2) 	# Feature selection rate
    TER <- round(100*sum(obj$yhat != yt) / nt, digits = 2)	# Test error rate
    res <- list(funCall = match.call(), class = "crda.demo", TER = TER, FSR = FSR, CT = CT )
    return(res)
}
