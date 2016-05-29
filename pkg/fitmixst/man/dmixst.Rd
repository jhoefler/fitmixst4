\name{dmixst}
\alias{dmixst}
\title{The value of the mixed densities of multivarite skew t distributions.}
\usage{
  dmixst(x, para)
}
\arguments{
  \item{x}{a possiblty multidimensional vector.}

  \item{para}{an object of the parameters of the mixture of
  multivariate skew t distributions. pro a vector for the
  mixture ratios, mu the values, sigma, delta, nu}
}
\value{
  numeric vector with the values of the mixture of the
  density functions.
}
\description{
  The value of the mixed densities of multivarite skew t
  distributions.
}
\examples{
x=c(1,2)
pro = 1; mu=1; Sigma=3; delta=2; nu=3;
para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
dmixst(x,para)
}
\keyword{density}
\keyword{function}

