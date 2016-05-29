

#' True Postive Count
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @return  true positive count
#' @description Area under the ROC curve. Calculated in various ways.
#' @author Andi Boeck
#' @export
TPC <- function(thres, y, y.hat)
{
  if(length(thres)==1) return(sum( (y.hat > thres) & (y == 1)))
  sapply(thres, TPC, y=y, y.hat=y.hat)
  
}



#' False Postive Count
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @export
#' 
FPC <- function(thres, y, y.hat)
{
  if(length(thres)==1) return(sum( (y.hat > thres) & (y == 0)))
  sapply(thres, FPC, y=y, y.hat=y.hat)
}

#' True Negative Count
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @export
#' 
TNC <- function(thres, y, y.hat)
{
  if(length(thres)==1) return(sum( (y.hat < thres) & (y == 0)))
  sapply(thres, TNC, y=y, y.hat=y.hat)
}

#' False Negative Count
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @export
#' 
FNC <- function(thres, y, y.hat)
{
  if(length(thres)==1) return(sum( (y.hat < thres) & (y == 1)))
  sapply(thres, FNC, y=y, y.hat=y.hat)
}

#' Sensitivity
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @export
#' 
sens <- function(thres, y, y.hat)
{
  TPC(thres=thres, y=y, y.hat=y.hat) / sum(y)
}

#' Specificity
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @param thres thershold where to split. Must be in range of \code{y.hat}
#' @export
#' 
spec <- function(thres, y, y.hat)
{
  1 - (FPC(thres=thres, y=y, y.hat=y.hat) / sum(y==0))
}


#' Number of concordant pairs
#' 
#' @param y.hat numeric. risk between 0 and 1 
#' @param y status yes=1, no=0 or dead=1, alive=0 
#' @description Concordant pairs
#' @export

Ncon <- function(y, y.hat) 
{
  #if(length(y)>1000) stop("y to big, use optimized code")
  a <- sum(sapply(y, function(x) (y < x)) * sapply(y.hat, function(x) (y.hat < x)))
  b <- sum(sapply(y, function(x) (y > x)) * sapply(y.hat, function(x) (y.hat >= x)))
  a + b
  # identisch zu
  #for(i in seq_along(y.hat)) {
  # for(j in seq_along(y.hat)) {
  #    N_c <- N_c + (y[i] < y[j] & y.hat[i] <  y.hat[j])
  #    N_c <- N_c + (y[i] > y[j] & y.hat[i] >=  y.hat[j])
  #  }
  #}
}


#' Number of relevant pairs
#'
#' @param y status yes=1, no=0 or dead=1, alive=0
#' @references Tutz, Analyse kategorialer Daten, p. 111ff 
#' @export
#' @description Relevant pairs
#' @examples
#'  outcome <- c(1,1,0,0); 
#'  Nrel(outcome)
Nrel <- function(y) 
{
  prod(table(y))*2
}

