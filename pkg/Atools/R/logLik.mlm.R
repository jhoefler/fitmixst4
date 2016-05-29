#' Log-likelihood of multivariate linear regression model
#' 
#' @method logLik mlm
#' @param object multivariate linear regression model fit with \code{\link{lm}}.
#'  of class \code{mlm}
#' @param ... not used.
#' @return log-lik at (unrestricted) maximum with df as attribute.
#' @author Andi Boeck
#' @import mvtnorm 
#' @export
#' @examples y <- cbind(rnorm(10), rnorm(10)); x <- 1:10;
#' mod <- lm(y~x)
#' logLik(mod)

logLik.mlm <- function(object,...)
{
  resids <- residuals(object)
  n <- nrow(resids)
  Sigma_ML <- crossprod(resids) /n
  ans <- sum(dmvnorm(resids, sigma=Sigma_ML, log=T))
  
  df <- length(coef(object)) + nrow(Sigma_ML) * (nrow(Sigma_ML) + 1) / 2
  attr(ans, "nobs") <- n
  attr(ans, "df") <- df
  class(ans) <- "logLik"
  ans
}

