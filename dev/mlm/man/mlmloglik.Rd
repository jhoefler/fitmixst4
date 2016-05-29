\name{mlmloglik}
\alias{mlmloglik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mlmloglik(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theta) 
{
    y <- mod$model[[1]]
    dim.y <- ncol(y)
    length.beta <- length(coef(mod))
    length.var <- dim.y
    length.corr <- (dim.y * (dim.y + 1))/2 - length.var
    betas <- matrix(theta[1:length.beta], ncol = dim.y, byrow = F)
    vars <- R2var(theta[(length.beta + 1):(length.beta + length.var)])
    corrs <- R2cor(theta[(length.beta + length.var + 1):(length.beta + 
        length.var + length.corr)])
    corr.mat <- diag(dim.y)
    corr.mat[lower.tri(corr.mat, diag = FALSE)] <- corrs
    corr.mat <- corr.mat + t(corr.mat)
    diag(corr.mat) <- 1
    var.diag <- diag(sqrt(vars))
    cov.mat <- var.diag \%*\% corr.mat \%*\% var.diag
    X <- model.matrix(mod)
    y.hat <- X \%*\% betas
    resids <- y - y.hat
    ans <- sum(mvtnorm::dmvnorm(resids, sigma = cov.mat, log = T))
    df <- length(coef(mod)) + dim.y * (dim.y + 1)/2
    attr(ans, "df") <- df
    -ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
