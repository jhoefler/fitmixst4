require(bbmle)

dmixnorm <- function(x, prop, mean, sd, log=TRUE)
{
    prop <- prop / sum(prop) 
    #in maply Beobachtungen n sind Spalten
    #Mischkomp in Zeilen    
    #dnorm(x, mean=mean, sd=sd)
    #do.call("dnorm", list(x=x, mean=mean, sd=sd))
    if(!log) 
    {
     return(colSums(t(mapply(dnorm, mean=mean, sd=sd, MoreArgs=list(x=x)))*prop))  
    }
    log(colSums(t(mapply(dnorm, mean=mean, sd=sd, MoreArgs=list(x=x)))*prop))
}

llmixnorm <- function(x, prop, mean, sd)
{
    
    sum(dmixnorm(x, prop,mean,sd))
}

llmixnorm1 <- function(args) 
{
    Ncomp <- length(args)/3
    
    prop <- args[1:Ncomp]
    prop <- scale.props(prop)
    mean <- args[(Ncomp+1):(2*Ncomp)]
    sd <- args[(2*Ncomp+1):(3*Ncomp)]
    -llmixnorm(y,prop,mean,sd)
}

llmixnorm2 <- function(strange.args) 
{
    Ncomp <- ceiling((length(strange.args))/3)
    args <- ll2args(strange.args)
    prop <- args[1:Ncomp]
    mean <- args[(Ncomp+1):(2*Ncomp)]
    sd <- args[(2*Ncomp+1):(3*Ncomp)]
    -llmixnorm(y,prop,mean,sd)
}

llmixnorm3 <- function(args) # mit log-odds
{
    Ncomp <- length(args)/3
    
    prop <- args[1:Ncomp]
    prop <- logodds2prop(prop)
    mean <- args[(Ncomp+1):(2*Ncomp)]
    sd <- args[(2*Ncomp+1):(3*Ncomp)]
    -llmixnorm(y,prop,mean,sd)
}


prop2logodds <- function(x)
{   # referenz ist erste kategroie
    x <- (x/x[1])[-1]
    log(x)
}
logodds2prop <- function(x)
{
    x <- c(1,exp(x))
    x/sum(x)
    
}

fitmixnorm <- function(y, start.method="kmeans", Ncomp = 2,...)
{
    
    km <- kmeans(y,Ncomp)
    means <- km$ce
    sds <- unlist(tapply(y, km$cl, sd))
    props <- unlist(table(km$cl) /length( km$cl))
    names(means) <- names(sds) <- names(props) <- NULL 
    
    args <- c(props, means, sds)
    strange.args <- args2fit(args)
    names(strange.args) <- strange.args2names(strange.args)
    names(args) <- args2names(args)
    
    #parnames(llmixnorm2) <- args2names(strange.args)
    #fit <- mle2(llmixnorm2, start=strange.args, data=list(y=y), vecpar=T,...)
    
    parnames(llmixnorm3) <- strange.args2names(strange.args)
    fit <- mle2(llmixnorm3, start=strange.args, data=list(y=y), vecpar=T,...)
    
    return(fit)
}


plot.mixnorm <- function(fit, x=NULL, ...)
{
    args <- ll2args(unlist(coef(fit)))
    Ncomp <- length(args)/3
    prop <- args[1:Ncomp]
    mean <- args[(Ncomp+1):(2*Ncomp)]
    sd <- args[(2*Ncomp+1):(3*Ncomp)]
    dat <- fit@data$y 
    if(is.null(x)) x <- seq(min(dat), max(dat), length=500)
    
    y.vals <- dmixnorm(x,prop,mean,sd,log=FALSE, ...)
    
    hist(dat, prob=T, breaks=length(dat)/10)
    abline(v=mean,lty=2, col="pink", lwd=3)
    lines(x, y.vals)
    
    adjust <- mean(sd) / (sd(dat)/sqrt(Ncomp))
    lines(density(dat, bw="SJ", adjust=adjust), col="red")
    
}

plot.mixnorm2 <- function(fit, x=NULL, ...)
{
    args <- (unlist(coef(fit)))
    Ncomp <- length(args)/3
    prop <- args[1:Ncomp]
    prop <- scale.props(prop)
    mean <- args[(Ncomp+1):(2*Ncomp)]
    sd <- args[(2*Ncomp+1):(3*Ncomp)]
    dat <- fit@data$y 
    if(is.null(x)) x <- seq(min(dat), max(dat), length=500)
    
    y.vals <- dmixnorm(x,prop,mean,sd,log=FALSE, ...)
    
    hist(dat, prob=T, breaks=length(dat)/10)
    abline(v=mean,lty=2, col="pink", lwd=3)
    lines(x, y.vals)
    
    adjust <- mean(sd) / (sd(dat)/sqrt(Ncomp))
    lines(density(dat, bw="SJ", adjust=adjust), col="red")
    
}

ll2args <- function(strange.args)
{
    
    Ncomp <- ceiling((length(strange.args))/3)
    
    
    strange.args.temp <- rev(strange.args)
    strange.args.temp[1:Ncomp] <- exp(strange.args.temp[1:Ncomp])
    strange.args <- rev(strange.args.temp)
    
    if(Ncomp==2) 
    {
        return(c(plogis(strange.args[1]), 1-plogis(strange.args[1]), strange.args[-1]))
    }
    prop <- double(Ncomp)
    prop[1] <- plogis(strange.args)[1]
    relfac <- plogis(strange.args[2:(Ncomp-1)])
    
    
    used.prop <- prop[1]
    for(i in seq_along(relfac)) 
    {
        prop[i+1] <- relfac[i] *(1-used.prop)
        used.prop <- used.prop + prop[i+1]
        prop[i+2] <- 1 - used.prop
    }
    args <- c(prop, strange.args[-c(1:(Ncomp-1))]) 
    return(args)
}

args2ll <- function(args)
{
    Ncomp <- floor(length(args)/3)
    args.temp <- rev(args)
    args.temp[1:Ncomp] <- log(args.temp[1:Ncomp])
    args <- rev(args.temp)
    
    if(Ncomp==2) 
    {
        return(c(qlogis(args[1]), args[-c(1,2)]))
    }
    prop <- args[1:Ncomp]
    relfac <- double(Ncomp-2)
    strange.prop <- double(Ncomp-1)
    
    strange.prop[1] <- qlogis(prop[1])
    used.prop <- 1- prop[1]
    for(i in 2:(Ncomp-1))
    {
        strange.prop[i] <- qlogis(prop[i] / used.prop)
        used.prop <- used.prop - prop[i]
        #strange.prop[i+1] <- qlogis(prop[i+1] / used.prop) 
    }
    strange.args <- c(strange.prop, args[-c(1:Ncomp)])
    return(strange.args)
}

strange.args2names <- function(strange.args)
{
    Ncomp <- ceiling((length(strange.args))/3)
    names <- c(paste("logOdd", 1:(Ncomp-1), sep=""), paste("mean", 1:Ncomp, sep=""), 
               paste("logSd", 1:Ncomp, sep=""))
    return(names)
}   

args2names <- function(args)
{
    Ncomp <- length(args)/3
    names <- c(paste("prop", 1:Ncomp, sep=""), paste("mean", 1:Ncomp, sep=""), 
               paste("sd", 1:Ncomp, sep=""))
    return(names)
} 

args2fit <- function(args)
{
    Ncomp <- floor(length(args)/3)
    args.temp <- rev(args)
    args.temp[1:Ncomp] <- log(args.temp[1:Ncomp])
    args <- rev(args.temp)
    
    strange.args <- c(prop2logodds(args[1:Ncomp]), args[-c(1:Ncomp)])
    return(strange.args)
}
