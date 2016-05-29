\name{spearman.ci}
\alias{spearman.ci}
\title{Asymptotic confidence intervals for correlation after Spearman}
\usage{
spearman.ci(x, y, alpha = 0.05)
}
\arguments{
  \item{x,y}{numeric vectors of same length}

  \item{alpha}{confidence level}
}
\value{
vector of length 3, cor, lower, upper
}
\description{
Asymptotic confidence intervals for correlation after
Spearman
}
\examples{
height <- c(1.56, 1.34, 1.45, 1.23, 1.30, 1.27)
weight <- c(23, 54, 34, 22, 20, 31)
spearman.ci(height, weight)
}
\author{
Andi Boeck
}
\references{
Bernard Rosner \emph{Fundamentals of Biostatistics}, 7th
Edition page 500. (ISBN-13: 9780538733496)
}

