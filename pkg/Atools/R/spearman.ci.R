#' Asymptotic confidence intervals for correlation after Spearman
#' 
#' @param x,y numeric vectors of same length
#' @param alpha confidence level
#' @return vector of length 3, cor, lower, upper
#' @author Andi Boeck
#' @references Bernard Rosner \emph{Fundamentals of Biostatistics}, 7th Edition
#' page 500. (ISBN-13: 9780538733496)
#' @export
#' @examples height <- c(1.56, 1.34, 1.45, 1.23, 1.30, 1.27)
#' weight <- c(23, 54, 34, 22, 20, 31)
#' spearman.ci(height, weight)


spearman.ci <- function(x,y, alpha=.05)
{
  
  r <- cor(x,y, method="spearman")
  n <- length(x)
  
  # rosner biostatbook, p.500
  #step 1
  p.hat <- rank(x) / (n+1)
  p.hat.star <- rank(y) / (n+1)
  
  
  h.star <- qnorm(p.hat)
  h.hat.star <- qnorm(p.hat.star)
  
  
  
  # step 2
  r.h <- cor(h.star, h.hat.star, method="pearson")
  
  #step 3
  r.cor.h <- r.h*(1+(1-r.h^2) / (2*(n-4)))
  
  
  #step 4
  fisher.trans <- function(x) 0.5*log( (1+x) / (1-x))
  fisher.back.trans <- function(x) (exp(2*x) -1) / (exp(2*x) +1)
  
  # step 5
  z.hat.h <- fisher.trans(r.cor.h)
  z.h <- z.hat.h +c(-1,1)*(qnorm(1-alpha/2)/sqrt(n-3))
  
  
  
  #step 6
  r.h <- fisher.back.trans(z.h)
  
  
  # step 7
  ans <- c(r, (6/pi)* asin(r.h/2))
  attr(ans, "names") <- c("rho", "lower", "upper")
  
  attr(ans, "Confidence level") <- 1-alpha
  attr(ans, "call") <- match.call()
  class(ans) <- "spearman.ci"
  
  # step 8
  # nothing to do...
  return(ans)
}

#' Print method for speraman.ci

#' @usage ##S3 method for class \code{spearman.ci}.
#' @param x Object of class \code{spearman.ci}.
#' @param digits Digits to show.
#' @param ... Not used.

print.spearman.ci <- 
  function(x, digits = max(3, getOption("digits") - 3),...)
{
  cat(paste("\nRank-correlation with asymptotic", 
            format(attr(x, "Confidence level")*100,digits=digits), 
            "% Confidence interval", "\n"))
  cat(format(x, digits = digits), "\n")
}





