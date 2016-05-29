\name{rmvn}
\alias{rmvn}
\title{Generates random observation for the normal distribution}
\usage{
  rmvn(n, p, mean = rep(0, p), cov = diag(p))
}
\arguments{
  \item{n}{number of random observations.}

  \item{p}{dimension of the skew t distribution.}

  \item{mean}{a mean vector with length p.}

  \item{cov}{a covariance matrix of dimension pxp.}
}
\value{
  gernerates random observations for the normal
  distribution.
}
\description{
  Generates random observation for the normal distribution
}
\examples{
mu=1; Sigma=1;
y <- rmvn(n=100,p=1,mean=mu,cov=Sigma)
}
\keyword{distribution.}
\keyword{observations}
\keyword{random}
\keyword{snormal}

