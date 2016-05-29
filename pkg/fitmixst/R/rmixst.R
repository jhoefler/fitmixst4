# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Gernates random observation for a mixed skew t distribution
#' @param n number of generated observations.
#' @param para the initial input parameters. Note: if you choose input parameters the kmeans method is disabled.
#' para a list of the input parameters. para = list (pro, mu, Sigma, delta, nu).
#' @return a vector of random observations.
#' @keywords generate random deviates.
#' @export
#' @examples
#' 
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)


rmixst<-function(n, para)
{
	pro <- para$pro
	mu <- para$mu
	Sigma <- para$Sigma
	delta <- para$delta
	nu <- para$nu
	
	g <- length(pro)
	y <- vector()
	
	if ( sum(pro) != 1)
	{
		stop("Sum of pro vector hast to be one")
	}
	for (i in 1:g)
	{
		y <- c(y,rmst(n*pro[i],p=1,mean=mu[i],cov=Sigma[i],nu=nu[i],del=delta[i]))
	}	
	return(y)
}

