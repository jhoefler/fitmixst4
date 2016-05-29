\name{rmst}
\alias{rmst}
\title{Generates random observation for the skew t distribution}
\usage{
  rmst(n, p, mean = rep(0, p), cov = diag(p), nu = 10,
    del = rep(0, p))
}
\arguments{
  \item{n}{number of random observations.}

  \item{p}{dimension of the skew t distribution.}

  \item{mean}{a mean vector with length p.}

  \item{cov}{a covariance matrix of dimension pxp.}

  \item{del}{a skew parameter of the length p.}

  \item{nu}{the degree of freedoms.}
}
\value{
  gernerates random observations for the skew t
  distribution.
}
\description{
  Generates random observation for the skew t distribution
}
\examples{
mu=1; Sigma=1; delta=3; nu=3;
y <- rmst(n=100,p=1,mean=mu,cov=Sigma,nu=nu,del=delta)
}
\keyword{distribution.}
\keyword{observations}
\keyword{random}
\keyword{skew}
\keyword{t}

