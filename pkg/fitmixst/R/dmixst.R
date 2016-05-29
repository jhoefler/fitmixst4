# TODO: Add comment
# 
# Author: jh
###############################################################################


#' The value of the mixed densities of multivarite skew t distributions.
#'
#' @param x a possiblty multidimensional vector.
#' @param para an object of the parameters of the mixture of multivariate skew t distributions.
#'pro a vector for the mixture ratios, mu the values, sigma, delta, nu
#' @return numeric vector with the values of the mixture of the density functions.
#' @keywords density function
#' @export
#' @examples
#' x=c(1,2)
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' dmixst(x,para)


dmixst <- function (x, para){
	pro<- para$pro
	mu <- para$mu
	Sigma <- para$Sigma
	delta <- para$delta
	nu <- para$nu
	.Call("file77c6a95fb7", x, pro, mu, Sigma, delta, nu, PACKAGE = "fitmixst")
}


