\name{auc}
\alias{auc}
\alias{aucHmisc}
\alias{aucPROC}
\alias{aucRLoop}
\alias{aucROCR}
\alias{aucRank}
\alias{aucRapply}
\alias{aucWilcox}
\title{AUC}
\usage{
auc(...)

aucPROC(y, y.hat, oriented = FALSE)

aucWilcox(y, y.hat, oriented = FALSE)

aucRank(y, y.hat, oriented = FALSE)

aucRLoop(y, y.hat, oriented = FALSE)

aucRapply(y, y.hat, oriented = FALSE)

aucROCR(y, y.hat, oriented = FALSE)

aucHmisc(y, y.hat, oriented = FALSE)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1.}

  \item{y}{numeric. Status yes=1, no=0 or dead=1, alive=0.}

  \item{oriented}{logical. If FALSE (default) returned AUC
  always \eqn{\ge} 0.5.}

  \item{...}{Pass arguments. If TRUE: Probability that
  higher y.hat correspond to higher y.}
}
\value{
numeric.
}
\description{
Area under the ROC curve. Calculated in various ways.
Fastest is \code{auc()}, currently set to \code{aucRank()}.
\code{RLoop()} only supports oriented=FALSE
}
\examples{
y.hat <- seq(0,1,length=100)
 y <- rbinom(100, size=1, prob=1-y.hat)
 aucRank(y=y, y.hat=y.hat, oriented=TRUE)
 aucRLoop(y=y, y.hat=y.hat, oriented=TRUE)
 aucPROC(y=y, y.hat=y.hat, oriented=TRUE)
 aucHmisc(y=y, y.hat=y.hat, oriented=TRUE)
 aucWilcox(y=y, y.hat=y.hat, oriented=TRUE)
 aucRapply(y=y, y.hat=y.hat, oriented=TRUE)
 aucPROC(y=y, y.hat=y.hat, oriented=FALSE)
 aucRank(y=y, y.hat=y.hat, oriented=FALSE)
 aucRLoop(y=y, y.hat=y.hat, oriented=FALSE)
 aucHmisc(y=y, y.hat=y.hat, oriented=FALSE)
 aucWilcox(y=y, y.hat=y.hat, oriented=FALSE)
 aucRapply(y=y, y.hat=y.hat, oriented=FALSE)
}
\author{
Andi Boeck
}

