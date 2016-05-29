\name{get.theta}
\alias{get.theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.theta(mod, REML = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
%%     ~~Describe \code{mod} here~~
}
  \item{REML}{
%%     ~~Describe \code{REML} here~~
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
function (mod, REML = FALSE) 
{
    resids <- residuals(mod)
    dim.y <- ncol(resids)
    n <- nrow(resids)
    nrow.beta <- nrow(coef(mod))
    Sigma_ML <- crossprod(resids)/n
    cor.mat <- cov2cor(Sigma_ML)
    betas <- c(coef(mod))
    if (rownames(coef(mod))[1] == "(Intercept)") {
        rownames(mod$coef)[1] <- "intercept"
    }
    names(betas) <- paste(rep(rownames(mod$coef), times = dim.y), 
        rep(1:dim.y, each = nrow.beta), sep = "")
    variances <- var2R(diag(Sigma_ML))
    names(variances) <- paste("sigma2.", 1:dim.y, sep = "")
    corrs <- cor2R(cor.mat[lower.tri(cor.mat, diag = FALSE)])
    names(corrs) <- paste("corr", 1:length(corrs), sep = "")
    ans <- as.list(c(betas, variances, corrs))
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
