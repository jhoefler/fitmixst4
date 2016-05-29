#' log-Likelihood of multivariate linear model
#' 
#' @param suitable theta parameter vector; at best taken form 
#' \link{get.theta}; corresponding model has to have name mod 
#' and be availible in parental frame (workspace).
#' @value Value of  negative log likelihood.
#' 

mlmloglik <-
function(theta)
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
