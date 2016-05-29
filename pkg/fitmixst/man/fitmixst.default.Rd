\name{fitmixst.default}
\alias{fitmixst.default}
\title{Fitted mixed model}
\usage{
  ##S3 method for class 'fitmixst' fitmixst(y,...)
}
\arguments{
  \item{y}{a multidimensional input vector.}

  \item{g}{the number of groups.}

  \item{rel.error}{the desired relative error. (default
  1e-5)}

  \item{itermax}{the maximum of iterations. (default 1000)}

  \item{method}{of finding initial values. (default
  "kmeans")}

  \item{para}{the initial input parameters. Note: if you
  choose input parameters the kmeans method is disabled.
  para a list of the input parameters. para = list (pro,
  mu, Sigma, delta, nu).}

  \item{...}{other paramters.}
}
\value{
  a object of the class fitmixst.
}
\description{
  Fitted mixed model
}
\examples{
pro = 1; mu=1; Sigma=3; delta=2; nu=3;
para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
y <- rmixst(100,para)
out <- fitmixst(y,g=2,method="kmeans")
}
\keyword{fit}
\keyword{mixed}
\keyword{skew}
\keyword{t.}

