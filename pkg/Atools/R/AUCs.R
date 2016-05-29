
#' AUC
#' 
#' @param y.hat numeric. risk between 0 and 1. 
#' @param y numeric. Status yes=1, no=0 or dead=1, alive=0.
#' @param oriented logical. If FALSE (default) returned AUC always \eqn{\ge} 0.5.
#' @param ... Pass arguments.
#' If TRUE: Probability that higher y.hat correspond to higher y.
#' @return numeric.
#' @description Area under the ROC curve. Calculated in various ways.
#' Fastest is \code{auc()}, currently set to \code{aucRank()}.
#' \code{RLoop()} only supports oriented=FALSE
#' @author Andi Boeck
#' @import pROC
#' @import ROCR
#' @export
#' @examples
#'  y.hat <- seq(0,1,length=100) 
#'  y <- rbinom(100, size=1, prob=1-y.hat) 
#'  aucRank(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucRLoop(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucPROC(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucHmisc(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucWilcox(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucRapply(y=y, y.hat=y.hat, oriented=TRUE)
#'  aucPROC(y=y, y.hat=y.hat, oriented=FALSE)
#'  aucRank(y=y, y.hat=y.hat, oriented=FALSE)
#'  aucRLoop(y=y, y.hat=y.hat, oriented=FALSE)
#'  aucHmisc(y=y, y.hat=y.hat, oriented=FALSE)
#'  aucWilcox(y=y, y.hat=y.hat, oriented=FALSE)
#'  aucRapply(y=y, y.hat=y.hat, oriented=FALSE)
#'  
auc <- function(...) aucRank(...)


#' @rdname auc
#' @export
aucPROC <- function(y, y.hat, oriented=FALSE) {
  if(oriented)
  {
    warning("only oriented=FALSE supported")
    return(NA)
  } 
  auc <- pROC:::auc(response=y, predictor=y.hat)[1]
  auc
}



#' @rdname auc
#' @export
aucWilcox <- function(y, y.hat, oriented=FALSE) {
  if(oriented)
  {
    warning("only oriented=FALSE supported")
    return(NA)
  } 
  n <- table(y)
  W <- wilcox.test(y.hat~y)$stat
  W / prod(n)
}  


#' @rdname auc
#' @export
aucRank <- function(y, y.hat, oriented=FALSE) {
  n <- table(y)
  W <- sum(rank(y.hat)[y==0]) - n[1]*(n[1]+1)/2
  auc <- 1 - W / prod(n)
  if(!oriented & ( auc < .5)) auc <- 1-auc
  auc
} 

#' @rdname auc
#' @export
aucRLoop <- function(y, y.hat, oriented=FALSE)
{
  
  N <- Nrel(y)
  
  N_c <- 0
  N_d <- 0
  
  for(i in seq_along(y.hat)) {
    for(j in seq_along(y.hat)) {
      
      
      # concordant pairs
      N_c <- N_c + (y[i] < y[j] & y.hat[i] <  y.hat[j])
      N_c <- N_c + (y[i] > y[j] & y.hat[i] >=  y.hat[j])
      
      # discordant pairs
      #N_d <- N_d + (y[i] < y[j] & y.hat[i] >  y.hat[j])
      #N_d <- N_d + (y[i] > y[j] & y.hat[i] <  y.hat[j])
    }
  }
  
  auc <- N_c / N
  if(!oriented & ( auc < .5)) auc <- 1-auc
  auc
}

#' @rdname auc
#' @export
aucRapply <- function(y, y.hat, oriented=FALSE) {
  
  N <- Nrel(y)
  N_c <- Ncon(y = y, y.hat = y.hat)
  
  auc <- N_c / N
  if(!oriented & ( auc < .5)) auc <- 1-auc
  auc
}

#' @rdname auc
#' @export
#' @import ROCR
aucROCR <- function(y, y.hat, oriented=FALSE) {
  auc <- ROCR::performance(ROCR::prediction(y.hat, y), "auc")@y.values[[1]]
  if(!oriented & ( auc < .5)) auc <- 1-auc
  auc
}



#' @rdname auc
#' @export
aucHmisc <- function(y,y.hat, oriented=FALSE)  {
  auc <- rcorr.cens(x=y.hat, S=y, outx=FALSE)[1]
  if(!oriented & ( auc < .5)) auc <- 1-auc
  auc
} 