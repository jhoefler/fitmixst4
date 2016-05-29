\name{logLik.mlm}
\alias{logLik.mlm}
\title{Log-likelihood of multivariate linear regression model}
\usage{
\method{logLik}{mlm}(object, ...)
}
\arguments{
  \item{object}{multivariate linear regression model fit
  with \code{\link{lm}}.  of class \code{mlm}}

  \item{...}{not used.}
}
\value{
log-lik at (unrestricted) maximum with df as attribute.
}
\description{
Log-likelihood of multivariate linear regression model
}
\examples{
y <- cbind(rnorm(10), rnorm(10)); x <- 1:10;
mod <- lm(y~x)
logLik(mod)
}
\author{
Andi Boeck
}

