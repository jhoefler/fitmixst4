# TODO: Add comment
# 
# Author: jh
###############################################################################


#' Internal fitting function
#'
#' @param y a multidimensional input vector.
#' @param g the number of groups.
#' @param error the desired relative error. (default 1e-5)
#' @param itermax the maximum of iterations. (default 1000)
#' @param pro1 the number of groups.
#' @param mu_neu the desired relative error. (default 1e-5)
#' @param Sigma_neu the maximum of iterations. (default 1000)
#' @param delta_neu of finding initial values. (default "kmeans")
#' @param nu_neu of finding initial values. (default "kmeans")
#' @param verbose of finding initial values. (default "kmeans")
#' @param ... other inputs
#' @return a object of the class fitmixst.
#' @keywords fit mixed skew t.
#' @export


mixstEstMCMC <- function(y, g, itermax, error, pro1, mu_neu, Sigma_neu, delta_neu, nu_neu, verbose=F,...)
{
	p <- 1
	n <- length(y)
	lik_ges <- lik_ges_neu <- 1	
	
	#function for estimateing nu
	nuf <- function(nu,pro2,e2,e1) log(nu/2)-digamma(nu/2)-1/sum(pro2)*sum((pro2*(e2-e1-1))) 
	
	iter <- 0
	#cl <- makeCluster(getOption("cl.cores", 4))
	
	while(iter < itermax && abs(lik_ges / lik_ges_neu-1) > error || iter == 0)
	{
		mu <- mu_neu
		delta_neu1<-delta <- delta_neu
		Sigma <- Sigma_neu
		nu <- nu_neu
		Delta <- delta
		
		Omega <- list()
		Omega_inv <- list()
		Lambda <- list()
		q <- d <- y_star <- e1 <- e2 <- e3 <- e4 <- t1 <- t2 <- t3 <- ew1 <- ew2 <- ew3 <- ew4 <- tmp <- S2 <-S3 <-pro2 <- matrix(NA,n,g)
		
		for(i in 1:g){
			#parSapply(1:g, function(i){	
			Omega[[i]]<-Sigma[[i]]+Delta[[i]]*t(Delta[[i]])
			Omega_inv[[i]] <- 1/Omega[[i]]
			Lambda[[i]]<-1-t(Delta[[i]])*Omega_inv[[i]]*Delta[[i]]
			pro2[,i]<-dmixst(y,list(pro=pro1[i],mu=mu[[i]],Sigma=Sigma[[i]],delta=delta[[i]],nu=nu[[i]]))/
					dmixst(y,list(pro=pro1,mu=mu,Sigma=Sigma,delta=delta,nu=nu))
			q[,i]<-(Delta[[i]])*c(Omega_inv[[i]])*(y-mu[[i]])
			d[,i]<-(y-mu[[i]])*c(Omega_inv[[i]])*(y-mu[[i]])
			y_star[,i]<-q[,i]*sqrt((nu[[i]]+p)/(nu[[i]]+d[,i]))
			
			# t1[,i]<-pt(q[,i]*sqrt((nu[[i]]+p+2)/(nu[[i]]+d[,i]))/sqrt(c(Lambda[[i]])),df=nu[[i]]+p+2)
			# t2[,i]<-pt(y_star[,i]/sqrt(c(Lambda[[i]])),df=nu[[i]]+p)
			
			# #calculating the expectations
			# ew2[,i] <- (nu[[i]]+p)/(nu[[i]]+d[,i])*t1[,i]/t2[,i]
			# ew1[,i] <- (ew2[,i] - log ((nu[[i]]+d[,i])/2) - (nu[[i]]+p)/(nu[[i]]+d[,i]) + digamma((nu[[i]]+p)/2))
			
			# tmp[,i]<-(nu[[i]]+p)/(nu[[i]]+d[,i])*pt(q[,i]/sqrt((nu[[i]]+d[,i])/(nu[[i]]+p+2)*
			# c(Lambda[[i]])),df=nu[[i]]+p+2)/pt(y_star[,i]/sqrt(c(Lambda[[i]])),df=nu[[i]]+p)
			
			# #moments of truncated t distribution
			# S2[,i]<-trunct.m1(q[,i],sqrt(((nu[[i]]+d[,i])/(nu[[i]]+p+2))*Lambda[[i]]),nu[[i]]+p+2)
			# S3[,i]<-trunct.m2(q[,i],sqrt(((nu[[i]]+d[,i])/(nu[[i]]+p+2))*Lambda[[i]]),nu[[i]]+p+2)
			
			# ew3[,i]<-tmp[,i]*S2[,i]
			# ew4[,i]<-tmp[,i]*S3[,i]
			
			a<-matrix(0,50,n)
			a<-sapply(1:n,function(j) rtmvt(50,mean=q[j,i],sigma=c(d[j,i]+nu[[i]])/(p+nu[[i]])*Lambda[[i]],
								df=round(nu[[i]]+p),lower=0,algorithm="rejection"))
			b<-matrix(0,50,n)
			for(k in 1:50)
				b[k,]<-sapply(1:n, function(j) (rgamma(1,shape=(nu[[i]]+2*p)/2,rate=((a[k,j]-q[j,i])*1/c(Lambda[[i]])*(a[k,j]-q[j,i])+d[j,i]+nu[[i]])/2)))
			
			#for(k in 1:50)
			#b[k,]=sapply(1:n, function(j) rgamma(1,(nu[[i]]+2*p)/2,((apply(a,2,mean))[j]-q[j,i])*1/c(Lambda[[i]])*((apply(a,2,mean))[j]-q[j,i])+d[j,i]+nu[[i]])/2)
			
			
			e1[,i]<-apply(log(b),2,mean)
			e2[,i]<-apply(b,2,mean)
			e3[,i]<-apply(a*b,2,mean)
			e4[,i]<-apply(a^2*b,2,mean)
			
			#updating parameters
			mu_neu[[i]]<-sum(pro2[,i]*(e2[,i]*y-Delta[[i]]*e3[,i]))/sum(pro2[,i]*e2[,i])
			
			delta_neu[[i]]<-sum(pro2[,i]*(y-mu_neu[[i]])*e3[,i])/sum(pro2[,i]*e4[,i])
			delta_neu1[[i]]<-sum(pro2[,i]*(y-mu_neu[[i]])*e3[,i])/sum(pro2[,i]*ew4[,i])
			
			Sigma_neu[[i]]<-1/sum(pro2[,i])*sum(pro2[,i]*(delta_neu[[i]]*e4[,i]*delta_neu[[i]]-
								(y-mu_neu[[i]])*e3[,i]*delta_neu[[i]]+(y-mu_neu[[i]])*(y-mu_neu[[i]])*e2[,i]-
								delta_neu[[i]]*e3[,i]*(y-mu_neu[[i]])))
			#Sigma_neu <- 1
			
			nu_neu[[i]]<-uniroot(function(nu) nuf(nu, pro2[,i], e2[,i], e1[,i]),interval=c(0.5,10e7))$root
			
			
			
		}
		
		#cat("\ne1 ", sum(abs(e1-ew1)/n), "e2 ", sum(abs(e2-ew2)/n), "e3 ", sum(abs(e3-ew3)/n), "e4 ", sum(abs(e4-ew4)/n))
		
		#cat(ew1)
		
		pro1<-apply(pro2,2,sum)/n
		
		iter <- iter+1
		
		lik_ges <- lik_ges_neu
		
		lik_ges_neu <- sum(log(dmixst(y,list(pro=pro1,mu=mu_neu,Sigma=Sigma_neu,delta=delta_neu,nu=nu_neu))))
		
		
		if(verbose) cat("Iteration: ", iter, "\nrelative error", abs(lik_ges / lik_ges_neu-1), "\nLoglikelihood: ", lik_ges_neu, 
					"\nmu_neu", mu_neu, "\tSigma_neu", Sigma_neu, "\tdelta", delta_neu, "\tnu", nu_neu, "\n\n")
		
	}
	
#sort by size for easy comparison of models
	l <- order(mu_neu, decreasing=T)
	mu_neu <- mu_neu[l]
	delta_neu <- delta_neu[l]
	Sigma_neu <- Sigma_neu[l]
	pro1 <- pro1[l]
	pro2 <- pro2[,l]
	
	list(coefficients = list(mu=mu_neu, Sigma=Sigma_neu, delta=delta_neu, nu=nu_neu, pro=pro1), iter=iter, logLik=lik_ges_neu, 
			posteriori=pro2)
	
}

