\name{mixstEstMCMC}
\alias{mixstEstMCMC}
\title{Internal fitting function}
\usage{
  mixstEstMCMC(y, g, itermax, error, pro1, mu_neu,
    Sigma_neu, delta_neu, nu_neu, verbose = F, ...)
}
\arguments{
  \item{y}{a multidimensional input vector.}

  \item{g}{the number of groups.}

  \item{error}{the desired relative error. (default 1e-5)}

  \item{itermax}{the maximum of iterations. (default 1000)}

  \item{pro1}{the number of groups.}

  \item{mu_neu}{the desired relative error. (default 1e-5)}

  \item{Sigma_neu}{the maximum of iterations. (default
  1000)}

  \item{delta_neu}{of finding initial values. (default
  "kmeans")}

  \item{nu_neu}{of finding initial values. (default
  "kmeans")}

  \item{verbose}{of finding initial values. (default
  "kmeans")}

  \item{...}{other inputs}
}
\value{
  a object of the class fitmixst.
}
\description{
  Internal fitting function
}
\keyword{fit}
\keyword{mixed}
\keyword{skew}
\keyword{t.}

