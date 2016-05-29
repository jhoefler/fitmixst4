\name{W}
\alias{W}
\title{Wilcoxon-statistic}
\usage{
W(y, y.hat)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}
}
\value{
numeric.
}
\description{
Simple wrapper around \code{wilcox.test}
}
\author{
Andi Boeck
}

