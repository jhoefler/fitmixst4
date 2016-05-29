# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Fitted mixed model
#' @usage ##S3 method for class 'fitmixst'
#' fitmixst(y,...) 
#' @param y a multidimensional input vector.
#' @param g the number of groups.
#' @param rel.error the desired relative error. (default 1e-5)
#' @param itermax the maximum of iterations. (default 1000)
#' @param method of finding initial values. (default "kmeans")
#' @param para the initial input parameters. Note: if you choose input parameters the kmeans method is disabled.
#' para a list of the input parameters. para = list (pro, mu, Sigma, delta, nu).
#' @param ... other paramters.
#' @return a object of the class fitmixst.
#' @keywords fit mixed skew t.
#' @examples
#' 
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)
#' out <- fitmixst(y,g=2,method="kmeans")

fitmixst.default <- function(y, g, rel.error = 1e-5, itermax=1000, method = "kmeans", para = NULL, ...)
{
	y <- as.vector(y)
	g <- as.integer(g)
	
	
	
	#using kmeans to get inital values
	if(method == "kmeans" || is.null(para)){
		
		init <- kmeans(y, g)
		
		pro1 <- init$size/length(y)
		mu_neu <- delta_neu <- Sigma_neu <- vector(length=g)
		for (j in 1:g) {
			mu_neu[j] <- init$centers[j, ]
			delta_neu[j] <- sign(sum((y[init$cluster == j] -(rep(mu_neu[j])))^3))
			Sigma_neu[j] <- sd(y[init$cluster == j])
			dimnames(Sigma_neu[j]) <- NULL
			names(mu_neu[j]) <- NULL
			names(delta_neu[j]) <- NULL
		}	
		nu_neu <-rep(4,g)
		
	}
	else{
		
		if(method =="self" && ! is.null(para)){
			pro1 <- para$pro
			mu_neu <- para$mu 
			Sigma_neu <- para$Sigma 
			delta_neu <- para$delta
			nu_neu <- para$nu
		}
		else stop("method and input parameters are not consistent")
		
	}
	
	est <- mixstEst(y=y, g=g, itermax=itermax, error=rel.error, pro=pro1, mu=mu_neu, Sigma=Sigma_neu, delta=delta_neu, nu=nu_neu)
	
  #cat(print(est))
	
	est$call <- match.call()
	
	class(est) <- "fitmixst"
  cat(print(est))
	est
	
}
