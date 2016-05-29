
data(swiss)
require(Atools)

mod <- lm(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~
        0+ Species, data=iris)


mod <- lm(cbind(Fertility, Agriculture, Infant.Mortality)~ 
    Examination + Education  + Catholic, data=swiss)

logLik(mod)



# bivariate log-lik
require(mvtnorm)
resids <- residuals(mod)
Sigma_ML <- crossprod(resids) / n
Sigma_wrong <- Sigma_ML
Sigma_wrong[c(2,3)] <- 0 # covariances to zero

# "correct" logLik
sum(dmvnorm(resids, sigma=Sigma_ML, log=T))


R2cor <- function(x)  2/(1+exp((-x)/5))-1
R2var <- exp
var2R <- log
cor2R <- function(x) -5*log( 2/(x+1) - 1)

get.theta <- function(mod, REML=FALSE) 
{
  resids <- residuals(mod)
  dim.y <- ncol(resids)
  n <- nrow(resids)
  nrow.beta <- nrow(coef(mod))
  Sigma_ML <- crossprod(resids)/n 
  cor.mat <- cov2cor(Sigma_ML)
  betas <- c(coef(mod))
  if(rownames(coef(mod))[1] == "(Intercept)")
  {
      rownames(mod$coef)[1] <- "intercept" 
  }
  
  names(betas) <- paste(rep(rownames(mod$coef), times=dim.y), rep(1:dim.y, each=nrow.beta), 
                        sep="")
  variances <- var2R(diag(Sigma_ML))
  names(variances) <- paste("sigma2.", 1:dim.y, sep="")
  corrs <- cor2R(cor.mat[lower.tri(cor.mat, diag=FALSE)]) 
  names(corrs) <- paste("corr", 1:length(corrs),sep="")
  ans <- as.list(c(betas, variances, corrs))
  ans
}

theta <- unlist(get.theta(mod))

mlmloglik <- function(theta)
{
  # theta: betas, variances, corrs
  y <- mod$model[[1]]
  dim.y <- ncol(y)
  length.beta <- length(coef(mod))
  length.var <- dim.y
  length.corr <- (dim.y * (dim.y+1)) / 2 - length.var

  betas <- matrix(theta[1:length.beta], ncol=dim.y, byrow=F)
  vars <- R2var(theta[(length.beta+1):(length.beta+length.var)])
  corrs <- R2cor(theta[(length.beta+length.var+1):
      (length.beta+length.var+length.corr)])
  corr.mat <- diag(dim.y)
  corr.mat[lower.tri(corr.mat, diag=FALSE)] <- corrs
  corr.mat <- corr.mat + t(corr.mat)
  diag(corr.mat) <- 1
  
  var.diag <- diag(sqrt(vars))
  cov.mat <- var.diag %*% corr.mat %*% var.diag
  
  X <- model.matrix(mod)
 
  y.hat <- X %*% betas
  resids <- y-y.hat
  ans <- sum(mvtnorm::dmvnorm(resids, sigma = cov.mat, log = T))
  df <- length(coef(mod)) + dim.y * (dim.y + 1)/2
  attr(ans, "df") <- df
  -ans
  
}

a <- mlmloglik(unlist(get.theta(mod)))




refit.mlm <- function(mod, start=NULL, ...)
{
  ####### define function with seperate arguments
  arglist <- get.theta(mod)
  args <- paste(names(arglist), collapse=",")
  string3 <- paste("mlmloglik(c(", args, "))", sep="")
  exp1 <- parse(text=string3)
  multipar <- function(x) a+b
  formals(multipar) <- arglist
  body(multipar) <- as.call(c(as.name("{"), exp1))  
  ######################################
  
  
  fit <- mle2(multipar, start=start,...)
  
  return(fit)
  
}

start <- get.theta(mod)


fit <- refit.mlm(mod, start=start)

pr <- profile(fit, which=13)

loglik2(1,2,3,4,5,6,7,8,9)
library(bbmle)

R2cor <- function(x)  2/(1+exp((-x)/5))-1
R2cor <- function(x)  x
R2var <- exp
R2var <- function(x) x
var2R <- log
var2R <- function(x) x
cor2R <- function(x) -5*log( 2/(x+1) - 1)
cor2R <- function(x) x

fit <- mle2(mlmloglik, start=unlist(get.theta(mod)), vecpar=T)

p1 <- profile(fit)
confint(fit, parm=13, method="uniroot")
plot(p1)


