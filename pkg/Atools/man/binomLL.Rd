\name{binomLL}
\alias{binomLL}
\title{Binomial log likelihood}
\usage{
binomLL(y, y.hat)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}
}
\value{
log likelihood, numeric.
}
\description{
little wrapper around \code{\link{dbinom}}
}
\author{
Andi Boeck
}

