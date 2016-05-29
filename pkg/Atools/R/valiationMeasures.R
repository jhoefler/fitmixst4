#' Cutoff with best Youden index
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return value of \code{y.hat} with best Youden index
#' @description Cutoff with best Youden index
#' @author Andi Boeck
#' @export
#' @examples
#'  y.hat <- seq(0,1,length=100) 
#'  y <- rbinom(100, size=1, prob=1-y.hat) 
#'  Ycutoff(y.hat = y.hat, y = y)

Ycutoff  <- function(y.hat, y) 
{
  res <- NA
  attr(res, "left=alive") <- NA
  if(length(unique(y)) == 1) return(res)
  n <- length(y)
  
  # direction left of thres = 0, right =1 ?
  direc <- (aucRank(y.hat=y.hat, y=y, oriented=TRUE) >= .5)
  names(direc) <- ""
  thres <- sort(unique(y.hat))
  if(length(thres)< 2) return(res)
  
  fun <- function(thres, y.hat, y) {
    # true positive count
    TPC <- TPC(y.hat = y.hat, y=y, thres = thres)
    
    # false negative count
    FNC <- FNC(y.hat = y.hat, y=y, thres = thres)
    
    # false positive count
    FPC <- FPC(y.hat = y.hat, y=y, thres = thres)
    
    # true negative count
    TNC <- TNC(y.hat = y.hat, y=y, thres = thres)
    
    Sens <- TPC / (TPC + FNC)
    Spec <- TNC / (TNC + FPC)
    
    return(Sens + Spec) # criterion to be maximized
  }
  
  Y <- sapply(thres, fun, y.hat=y.hat, y=y)
  #cat(table(thres), "\n")
  
  res <- thres[which.max(Y)]
  if(!direc) res <- thres[which.min(Y)]
  attr(res, "left=alive") <- direc
  
  return(res)
}


#' Binomial log likelihood
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return log likelihood, numeric.
#' @description little wrapper around \code{\link{dbinom}}
#' @author Andi Boeck
#' @export
binomLL <- function(y, y.hat) sum( dbinom(x=y, size=1, prob=y.hat, log=T))




#' Nagelkerke R-square measure
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return R-square, standardized version.
#' @description likelihood based R-square
#' @author Andi Boeck
#' @export

NagR2 <- function(y, y.hat) {
  # identisch zu R2 in rms::val.prob 
  # mit dbinom schneller als Likelihood ausgeschrieben
  pi0 <- mean(y)
  n <- length(y)
  
  L0 <- binomLL(y=y, y.hat = pi0)
  L <- binomLL(y=y, y.hat=y.hat)
  
  
  maxR2 <- 1 - exp(L0*(2/n))
  R2 <- 1- exp((L0 - L)*(2/n))
  R2nag <- R2 / maxR2
  
  R2nag
}







#' Bierer score / MSPE for non-cencored outcomes
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param scaled If TRUE scaled version is used. 
#' @return Brier score numeric.
#' @description with no survival outcome, just the mean squared predicted 
#' error (MSPE). 
#' @details The scaling depends 
#' on the predictions \code{y.hat} and not only
#' on the actual outcome \code{y} as in \code{NagR2}. This limits the use 
#' of the scaled version to assess different
#' models on external data.
#' @author Andi Boeck
#' @export
Brier <- function(y, y.hat, scaled=FALSE) { 
  B <- mean((y-y.hat)^2)
  if(scaled)
  {
    mean.y.hat <- mean(y.hat)
    Bmax <- mean.y.hat * (1 - mean.y.hat)
    B <- 1 - B/Bmax
  }
  return(B)
}






#' Discrimination slope
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return postive, numeric.
#' @description Difference of the mean predictions for the different outcomes.
#' @details If positve result: average (mean) prediction of y=1 group is higher.
#' @author Andi Boeck
#' @export
dSlope <- function(y, y.hat) {
  # discrimination slope = diff in means steyerberg 2012 Epidemiology table 1
  mean(y.hat[y==1]) - mean(y.hat[y==0])
}






#' Calibration slope
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Calibration slope, 
#' @author Andi Boeck
#' @export
cSlope <- function(y, y.hat) {
  # slope of regreesoin y.hat
  coef(glm(y ~ qlogis(y.hat), family=binomial))[2]
}






#' Calibration in-the-large
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Calibration in-the-large
#' @details Mean prediction minus mean outcome.
#' @author Andi Boeck
#' @export
cLarge <-  function(y, y.hat) {
  # calibration in the large teyerberg 2012 Epidemiology table 1
  mean(y.hat) - mean(y)
}









#' Hosmer-Lemeshow statistic
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Groups by deciles.
#' @author Andi Boeck
#' @export
#' @import Hmisc
#' @details Chi-Square typ measure.
Hosmer <- function(y, y.hat) {
  N <- length(y.hat)
  mybreaks <- quantile(y.hat, seq(0,1,by=.1))
  grouping <- cut2(y.hat, cuts=mybreaks)
  n <- table(grouping)
  observed <- tapply(y, grouping, sum)
  expected <- tapply(y.hat, grouping, sum)
  mypi <- observed / n
  
  
  stat <- sum((observed - expected)^2 / (n*(mypi)*(1-mypi)) )
  return(stat)
}






#' Hosmer-Lemeshow statistic
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param g how many groups.
#' @return numeric.
#' @description Differt version to \code{Hosmer}.
#' @author Andi Boeck
#' @export
#' @import Hmisc
#' @details Chi-Square typ measure.
Hosmer2 <- function (y, y.hat, g = 10)
{
  cutyhat <- cut(y.hat, breaks = quantile(y.hat, probs = seq(0,
                                                             1, 1/g)), include.lowest = T)
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - y.hat, y.hat) ~ cutyhat)
  chisq <- sum((obs - expect)^2/expect)
  P <- 1 - pchisq(chisq, g - 2)
  c("X^2" = chisq)#, Df = g - 2, "P(>Chi)" = P)
}





#' Wilcoxon statistic by hand
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Wilcoxon statistic by hand.
#' @author Andi Boeck
#' @export
#' @details Just to know how R does calculate it.
WilcoxHand <- function(y, y.hat)
{
  # ?Wilcoxon # RHilfe
  N <- 0
  x <- y.hat[y==0]
  y2 <- y.hat[y==1]
  for(i in seq_along(x)) {
    for(j in seq_along(y2)) {
      
      N <- N + (x[i] <= y2[j])
      
      
    }
  }
  N
}



#' t-statistic
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Simple wrapper around \code{t.test}
#' @author Andi Boeck
#' @export
tstat <- function(y , y.hat) t.test(y.hat ~ y)$statistic


#' Wilcoxon-statistic
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Simple wrapper around \code{wilcox.test}
#' @author Andi Boeck
#' @export
W <- function(y, y.hat) wilcox.test(y.hat ~ y)$stat




#' Mean of sqaured Pearson-residuals
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Pearson Residual
#' @author Andi Boeck
#' @export
meanPearson <- function(y, y.hat) {
  mean ( ((y - y.hat)^2) / (y.hat*(1-y.hat)))
}




#' Calibration intercept
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return numeric.
#' @description Intercept of regreestion with slope fixed to one, 
#' measures calibration-in-the-large
#' @author Andi Boeck
#' @export
cInt <- function(y, y.hat) {
  coef(glm(y ~ 1 + offset(qlogis(y.hat)), family=binomial))[1]
}




#' Calculate netbenefit curves of risk predictions
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param Ngrid Number of grid points. Defaults to NULL where all distint \code{y.hat}
#'  are used.
#' @return data.frame with columns netben, all and threshold
#' @author Andi Boeck
#' @export
#' @examples
#' require(ggplot2)
#' set.seed(123)
#' risk <- runif(30, .1, .9);
#' outcome <- rbinom(30, size=1, prob=risk); 
#' bla <- netben(risk, outcome)
#' 
#' bla2 <- ggplot(bla, aes(x=threshold*100)) +
#'  geom_hline(aes(yintercept=0, linetype="No biopsies")) +
#'  geom_line(aes(y=all, linetype="Biopsy all"),alpha=.7, size=.4) +
#'  geom_line(aes(y=netben, linetype="Biopsy with risk threshold"),alpha=.7, size=.4) +
#'  xlab("Risk (%)") + ylab("Net benefit") +
#'  scale_linetype_discrete(guide = guide_legend(title=""))
#' print(bla2)

netben <- function(y.hat, y, Ngrid=NULL) 
{
  
  n <- length(y.hat)  
  thresholds <- c(0,sort(unique(y.hat)),1)
  if(!is.null(Ngrid)) 
  {
    thresholds <- c(0,seq(min(y.hat), max(y.hat), length=Ngrid),1)
  }
  
  fun <- function(thres, y.hat, y) {
    # true positive count
    TPC <- TPC(y.hat = y.hat, y=y, thres = thres)
    
    
    # false positive count
    FPC <- n - TPC
    
    cterm <- thres / (1 - thres)
    netben <- TPC/n - FPC/n * cterm
    
    # all to biopsy   %y-%Non-y × c/(1-c).
    all <- mean(y==1) - mean(y==0) * cterm
    return(data.frame(netben=netben, all=all, threshold=thres))
  }
  res <- as.data.frame(t(sapply(thresholds, fun, y.hat=y.hat, y=y)))
  data.frame(lapply(res[-nrow(res),], unlist))
}


#' Combine different validation measures in \code{data.frame}.
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return data.frame with several measures as columns
#' @author Andi Boeck
#' @import plyr
#' @export
#' @examples
#' outcome <- c(1,1,0,0); marker <- c(.3, .2, .1,.8);
#' oisse(marker, outcome)
oisse <- function(y, y.hat)
{
  data.frame(auc = aucRank(y = y, y.hat = y.hat),
             Brier = Brier(y = y, y.hat = y.hat),
             #binomLL = binomLL(y = y, y.hat = y.hat),
             #cLarge = cLarge(y = y, y.hat = y.hat),
             #cSlope = cSlope(y = y, y.hat = y.hat),
             #dSlope = dSlope(y = y, y.hat = y.hat),
            # Homser = Hosmer(y = y, y.hat = y.hat),
             NagR2 = NagR2(y = y, y.hat = y.hat))
}
             
#' Prepare \code{data.frame} for ROC-curve
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @return \code{data.frame}
#' @description Cutoff with best Youden index
#' @author Andi Boeck
#' @export
#' @details If low values of predictor/risk/\code{y.hat} correspond to event of
#'  interest coded with \code{y=1}, this will result in an ROC-curve below
#'  the diagonal. 
#' @examples
#'  n <- 1e2
#'  y.hat <- seq(0, 1, length=n) 
#'  y <- rbinom(n, size = 1, prob = y.hat) 
#'  rocdat <- ROCdat(y.hat = y.hat, y = y)
#'  
#'  with(rocdat, plot(x= 1 - Specificity, y = Sensitivity, type = "l"))
#'  (a <- Ycutoff(y.hat=y.hat, y))
#'  (p <- rocdat[which(rocdat$thres == a),])
#'  points(x=1-p$Spec, y=p$Sens, col="red")
#'  
#'  # reverse direction
#'  y <- rbinom(n, size = 1, prob = 1- y.hat) 
#'  rocdat <- ROCdat(y.hat = y.hat, y = y)
#'  with(rocdat, plot(x = 1 - Specificity, y = Sensitivity, type = "l"))

ROCdat  <- function(y.hat, y) 
{
  rev(thres <- c(-Inf,sort(unique(y.hat)),Inf))
  res <- data.frame(Sensitivity = sens(thres=thres, y.hat=y.hat, y=y),
                    Specificity = spec(thres=thres, y.hat=y.hat, y=y),
                    thres = thres)
  res
}

