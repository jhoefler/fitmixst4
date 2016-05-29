require(lme4)
sleepstudy$bla <- rnorm(nrow(sleepstudy))
sleepstudy$bla2 <- gl(5, nrow(sleepstudy)/5)
(fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE))

(fm1 <- lmer(Reaction ~ Days + bla + bla2 +(1|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + bla+ (1|Subject), sleepstudy))

mixed(Reaction ~ Days + bla + bla2+ (1|Subject), data = sleepstudy)
anova(fm1,fm2)


theta <- c(fixef(fm1), c(unlist(ranef(fm1))), log(c(VarCorr(fm1)[[1]])),
    log(summary(fm1)@sigma^2))

X <- fm1@X
U <- t(as.matrix(fm1@Zt))
y <- fm1@y
A <- diag(nrow(U))


ll <- function(theta,X,U,y,A)
{ 
  # notation nach Fahrmeir Regressionsbuch
  length.beta <- ncol(X)
  length.gamma <- ncol(U)
  
  beta <- theta[(1:length.beta)]
  gamma <- theta[((length.beta+1):(length.beta+length.gamma))]
  
  sigma.b <- theta[(length.beta+length.gamma+1)] 
  sigma <- theta[(length.beta+length.gamma+2)]
  
  n <- nrow(U)
  
  
  R <- diag(n)
  G <- t(U) %*% (sigma.b*A) %*% (U)
  
  ll_fixed <-
  
  ll <- -1/2 * crossprod((y - X%*%beta - U%*%gamma)) - # crosp wenn R diag(n)
    -1/2 * t(gamma) %*% solve(G) %*% gamma
  
  return(ll)
  
#   # evaluate covariance matrix for y
#   V <- U%*%A%*%t(U)*sigma.bˆ2 + diag(n)*sigmaˆ2
#   L <- chol(V) # L’L=V
#   # transform dependent linear model to indep.
#   y <- backsolve(L,y,transpose=TRUE)
#   X <- backsolve(L,X,transpose=TRUE)
#   b <- coef(lm(y˜X-1)) # estimate fixed effects
#   # evaluate log likelihood
#   logLik <- -n/2*log(2*pi) -
#     sum (log(diag(L))) -
#     sum((y-X%*%b)ˆ2)/2
#   attr(logLik,"fixed") <- b # allow retrieval of beta
#   logLik
}


ll(theta,X,U,y,A)
theta2 <- theta
theta2[21] <- .001
ll(theta2,X,U,y,A)

require(mgcv)
G <- t(U) %*% A %*% (U)
mod <- gam(y ~ 0 + X + U, paraPen=list(U=list(G)), method="REML")

d <- profile(fm1)

llmm <- function(theta,X,U,y,A,smin)
{ 
  # notation nach Regressionbuch Kasten S. 260
  length.beta <- ncol(X)
  length.gamma <- ncol(U)
  
  betas <- theta[(1:length.beta)]
  gammas <- theta[((length.beta+1):(length.beta+length.gamma))]
  
  sigma.b <- max(exp(theta[(length.beta+length.gamma+1)]), smin)
  sigma <- max(exp(theta[(length.beta+length.gamma+2)]),smin)
  
  N <- nrow(U)
  
  pred.random <- U %*% gammas
  
  fitted <- X %*% betas + pred.random
  eps <- y - fitted
  
  dat <- cbind(eps, pred.random)
  
  Sigma <- diag(c(sigma, sigma.b))
  
  ans <- sum(dmvnorm(dat, sigma=Sigma), log=T)
  attr(ans, "fixed") <- betas
  attr(ans, "random") <- gammas
  attr(ans, "sigma") <- (sigma)
  attr(ans, "sigma.b") <- (sigma.b)
  
  ans
}

llmm(theta,X,U,y,A,smin=.00001)
theta2 <- theta
theta2[21] <- .1
llmm(theta2,X,U,y,A)
start <-  rnorm(22, theta, sd=2)
llmm(start,X,U,y,A,smin=.00001)

start[1] <- 100
theta[1] <- 100
res <- optim(theta, llmm, X=X, U=U, y=y, A=A,smin=0, method="CG",
             control=list(fnscale=-1, maxit=10000))
llmm(res$par,X,U,y,A,smin=.0001)
plot(attr(llmm(res$par,X,U,y,A,smin=.0001), "random"), ranef(fm1)[[1]][,1])
abline(c(0,1))

require(synbreed)
data(maize)
maizeC <- codeGeno(maize)

# pedigree-based (expected) kinship matrix
K <- kin(maizeC,ret="kin",DH=maize$covar$DH)

########################################
# notation nach dem doktorarbeits-schnipsel
# braucht
require(lme4)
require(mvtnorm)
dat <- subset(sleepstudy, Subject %in% c(330:335))
(fm1 <- lmer(Reaction ~ Days + (1|Subject), dat, REML=FALSE))
X <- fm1@X
Z <- t(as.matrix(fm1@Zt))
y <- fm1@y

LL <-function(beta1=fixef(fm1)[1], beta2=fixef(fm1)[2], g1=fm1@ranef[1], 
              g2=fm1@ranef[2], g3=fm1@ranef[3],
              g4=fm1@ranef[4], g5=fm1@ranef[5], g6=fm1@ranef[6], 
              sigmaE=6.945435, sigmaB=log(VarCorr(fm1)[[1]][1]))
{
 gamma <- c(g1, g2, g3, g4, g5, g6)
 beta <- c(beta1, beta2)
 sigmaE <- exp(sigmaE); sigmaB <- exp(sigmaB)
 G <- diag(length(gamma))*sigmaB
 R <- diag(length(y))*sigmaE
 
 resp <- c(gamma, y - X%*%beta - Z%*%gamma)
 bigCov <- as.matrix(bdiag(G, R))
 #bigCov[bigCov>100] <- 100
 cat(diag(bigCov), "\n")
 
 -sum(dmvnorm(resp, sigma=bigCov, log=T))
}
 
LL()
require(bbmle)
fit <- mle2(LL, control=list(trace=T, REPORT=50))



length.beta <- ncol(X)
length.gamma <- ncol(U)

beta <- theta[(1:length.beta)]
gamma <- theta[((length.beta+1):(length.beta+length.gamma))]

sigma.b <- theta[(length.beta+length.gamma+1)] 
sigma <- theta[(length.beta+length.gamma+2)]

n <- nrow(U)


R <- diag(n)
G <- t(U) %*% (sigma.b*A) %*% (U)

