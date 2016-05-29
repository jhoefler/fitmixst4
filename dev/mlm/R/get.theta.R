#' get.theta
#' 
#' @param mod model of class \code{mlm}.
#' @param REML not effect til now, only ML availible
#' @return named list with ML-parameters

get.theta <-
function(mod, REML=FALSE) 
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
