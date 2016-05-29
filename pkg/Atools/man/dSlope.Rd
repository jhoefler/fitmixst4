\name{dSlope}
\alias{dSlope}
\title{Discrimination slope}
\usage{
dSlope(y, y.hat)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}
}
\value{
postive, numeric.
}
\description{
Difference of the mean predictions for the different
outcomes.
}
\details{
If positve result: average (mean) prediction of y=1 group
is higher.
}
\author{
Andi Boeck
}

