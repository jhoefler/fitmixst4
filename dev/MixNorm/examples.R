source("functions.R") # working directory here !

### example 1
require(flexmix)
set.seed(1)
y <- c(rnorm(200, 0,1), rnorm(200, -5, 1), rnorm(200, 5,1))


fit1 <- fitmixnorm(y,Ncomp=4, method="L-BFGS-B")
fit2 <- fitmixnorm(y,Ncomp=3, method="L-BFGS-B")
fit3 <- flexmix(y~1,cluster=kmeans(y,3)$cl)
fit4 <- flexmix(y~1,cluster=apply(rmultinom(length(y),1, runif(3)),2,function(x) which(x==1)))
logLik(fit4)

anova(fit1, fit2)
BIC(fit1, fit2)

plot.mixnorm2(fit1)
plot.mixnorm2(fit2)

summary(fit2)
pr <- profile(fit2, try_harder=TRUE)
(cis <- apply(confint(pr),2,ll2args)) # ci on understandable parameter scale


#############################################################################
### example 2
set.seed(1)
y <- c(rnorm(800, -5,1), rnorm(400, -3, 2), rnorm(200, 0,3),
       rnorm(100, 2, 2), rnorm(500, 5, .5))


fit <- fitmixnorm(y,Ncomp=5, method="L-BFGS-B", control=list(maxit=10e6,
        REPORT=50, trace=6)) # werden 465 iterations
fit@details
summary(fit)
ll2args(coef(fit))

plot.mixnorm(fit)


fit2 <- fitmixnorm(y,Ncomp=5, method="BFGS", control=list(maxit=10e6,
                    REPORT=50, trace=6))
ll2args(coef(fit2))
plot.mixnorm(fit2)

fit3 <- fitmixnorm(y,Ncomp=3, method="L-BFGS-B", control=list(maxit=10e6,
                             REPORT=50, trace=6))
ll2args(coef(fit3))
plot.mixnorm(fit3)

BIC(fit, fit2, fit3)

