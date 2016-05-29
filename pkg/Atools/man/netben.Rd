\name{netben}
\alias{netben}
\title{Calculate netbenefit curves of risk predictions}
\usage{
netben(y.hat, y, Ngrid = NULL)
}
\arguments{
  \item{y.hat}{numeric. risk between 0 and 1}

  \item{y}{status yes=1, no=0 or dead=1, alive=0}

  \item{Ngrid}{Number of grid points. Defaults to NULL
  where all distint \code{y.hat} are used.}
}
\value{
data.frame with columns netben, all and threshold
}
\description{
Calculate netbenefit curves of risk predictions
}
\examples{
require(ggplot2)
set.seed(123)
risk <- runif(30, .1, .9);
outcome <- rbinom(30, size=1, prob=risk);
bla <- netben(risk, outcome)

bla2 <- ggplot(bla, aes(x=threshold*100)) +
 geom_hline(aes(yintercept=0, linetype="No biopsies")) +
 geom_line(aes(y=all, linetype="Biopsy all"),alpha=.7, size=.4) +
 geom_line(aes(y=netben, linetype="Biopsy with risk threshold"),alpha=.7, size=.4) +
 xlab("Risk (\%)") + ylab("Net benefit") +
 scale_linetype_discrete(guide = guide_legend(title=""))
print(bla2)
}
\author{
Andi Boeck
}

