\name{Ycutoff}
\alias{Ycutoff}
\title{Cutoff with best Youden index}
\usage{
Ycutoff(y.hat, y)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}
}
\value{
value of \code{y.hat} with best Youden index
}
\description{
Cutoff with best Youden index
}
\examples{
y.hat <- seq(0,1,length=100)
 y <- rbinom(100, size=1, prob=1-y.hat)
 Ycutoff(y.hat = y.hat, y = y)
}
\author{
Andi Boeck
}

