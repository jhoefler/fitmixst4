% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/mmixstEst2.R
\name{mmixstEst2}
\alias{mmixstEst2}
\title{Internal fitting function}
\usage{
mmixstEst2(y, g, itermax, error, pro1, mu_neu, Sigma_neu, delta_neu, nu_neu,
  verbose = F, mcmc = F, ...)
}
\arguments{
\item{y}{a multidimensional input vector.}

\item{g}{the number of groups.}

\item{itermax}{the maximum of iterations. (default 1000)}

\item{error}{the desired relative error. (default 1e-5)}

\item{pro1}{the number of groups.}

\item{mu_neu}{the desired relative error. (default 1e-5)}

\item{Sigma_neu}{the maximum of iterations. (default 1000)}

\item{delta_neu}{of finding initial values. (default "kmeans")}

\item{nu_neu}{of finding initial values. (default "kmeans")}

\item{verbose}{of finding initial values. (default "kmeans")}

\item{mcmc}{calculate moments of truncated t distriubtion using mcmc algoirthm.}

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

