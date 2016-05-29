\name{Brier}
\alias{Brier}
\title{Bierer score / MSPE for non-cencored outcomes}
\usage{
Brier(y, y.hat, scaled = FALSE)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}

  \item{scaled}{If TRUE scaled version is used.}
}
\value{
Brier score numeric.
}
\description{
with no survival outcome, just the mean squared predicted
error (MSPE).
}
\details{
The scaling depends on the predictions \code{y.hat} and not
only on the actual outcome \code{y} as in \code{NagR2}.
This limits the use of the scaled version to assess
different models on external data.
}
\author{
Andi Boeck
}

