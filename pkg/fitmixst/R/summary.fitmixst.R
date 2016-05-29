# TODO: Add comment
# 
# Author: jh
###############################################################################


#' print for fitmixst
#'
#' @usage ##S3 method for class 'fitmixst'
#' summary(object,...)
#' @param object of the fitmixst class.
#' @return summary of the object
#' @keywords summary of the object
#' @export
#'

summary.fitmixst <- function (object, ...)
{
	cat("Call :\n")
	print(object$call)
	cat("\nIterations: \n")
	print(object$iter)
	cat("\nLogLikelihood: \n")
	print(object$logLik)
}

