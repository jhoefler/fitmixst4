# TODO: Add comment
# 
# Author: jh
###############################################################################


#' print for fitmixst
#'
#' @usage ##S3 method for class 'fitmixst'
#' print(object,...)
#' @param object of the fitmixst class.
#' @return the loglikhood of the object
#' @keywords printed object
#' @export
#' 
print.fitmixst <- function (x, ...)
{
	cat("Call :\n")
	print(x$call)
	cat("\nCoefficients: \n")
	print(x$coefficients)
	cat("\nIterations: \n")
	print(x$iter)
	cat("\nLogLikelihood: \n")
	print(x$logLik)
}

