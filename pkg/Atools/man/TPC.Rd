\name{TPC}
\alias{TPC}
\title{True Postive Count}
\usage{
TPC(thres, y, y.hat)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}

  \item{thres}{thershold where to split. Must be in range
  of \code{y.hat}}
}
\value{
true positive count
}
\description{
Area under the ROC curve. Calculated in various ways.
}
\author{
Andi Boeck
}

