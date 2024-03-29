\name{rmixst}
\alias{rmixst}
\title{Gernates random observation for a mixed skew t distribution}
\usage{
  rmixst(n, para)
}
\arguments{
  \item{n}{number of generated observations.}

  \item{para}{the initial input parameters. Note: if you
  choose input parameters the kmeans method is disabled.
  para a list of the input parameters. para = list (pro,
  mu, Sigma, delta, nu).}
}
\value{
  a vector of random observations.
}
\description{
  Gernates random observation for a mixed skew t
  distribution
}
\examples{
pro = 1; mu=1; Sigma=3; delta=2; nu=3;
para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
y <- rmixst(100,para)
}
\keyword{deviates.}
\keyword{generate}
\keyword{random}

