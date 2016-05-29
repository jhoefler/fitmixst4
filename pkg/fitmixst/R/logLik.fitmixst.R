# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Loglikelihood for fitmixst
#'
#' @usage ##S3 method for class 'fitmixst'
#' logLik(object,...)
#' @param object of the fitmixst class.
#' @return the loglikhood of the object
#' @keywords log likelihood
#' @export


logLik.fitmixst <- function (object, ...)
{
	cat("LogLikelihood :\n")
	object$logLik
}

