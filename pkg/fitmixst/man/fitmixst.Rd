\name{fitmixst}
\alias{fitmixst}
\title{Fitted mixed model}
\usage{
  fitmixst(y, ...)
}
\arguments{
  \item{y}{a multidimensional input vector.}

  \item{...}{the rest}
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

