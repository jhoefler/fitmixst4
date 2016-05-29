# TODO: Add comment
# 
# Author: jh
###############################################################################


mixstEstMCMC <- function(y, g, itermax, error, pro1, mu_neu, Sigma_neu, delta_neu, nu_neu, verbose=F,...)
{
	p <- dim(mean)[2]
	n <- length(y)
	lik_ges <- lik_ges_neu <- 1	
	
	#function for estimateing nu
	nuf <- function(nu,pro2,e2,e1) log(nu/2)-digamma(nu/2)-1/sum(pro2)*sum((pro2*(e2-e1-1))) 
	
	iter <- 0
	
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
		q <- array(NA,c(n,p,g))
		
		for(i in 1:g){
			
			Delta=diag(delta[[i]])
			Omega[[i]]<-Sigma[[i]]+(Delta[[i]])%*%t(Delta[[i]])
			Omega_inv[[i]] <- solve(Omega[[i]])
			Lambda[[i]]<-1-t(Delta[[i]])%*%Omega_inv[[i]]%*%Delta[[i]]
			
			for(j in 1:n){
				
				
				# TODO: multivariate density
				#pro2[,i]<-dmixst(y,list(pro=pro1[i],mu=mu[[i]],Sigma=Sigma[[i]],delta=delta[[i]],nu=nu[[i]]))/
				#		dmixst(y,list(pro=pro1,mu=mu,Sigma=Sigma,delta=delta,nu=nu))
				# q is an nxp matrix
				q[j,i]<-(Delta[[i]])%*%(Omega_inv[[i]])%*%(y[j,]-mu[[i]])
				d[j,i]<-(y[j,]-mu[[i]])%*%(Omega_inv[[i]])%*%(y[j,]-mu[[i]])
				y_star[j,i]<-q[j,i]*sqrt((nu[[i]]+p)/(nu[[i]]+d[j,i]))
				
				
				
				a<-matrix(0,50,n)
				a<-sapply(1:n,function(j) rtmvt(50,mean=q[j,i],sigma=c(d[j,i]+nu[[i]])/(p+nu[[i]])*Lambda[[i]],
									df=round(nu[[i]]+p),lower=0,algorithm="gibbs"))
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


