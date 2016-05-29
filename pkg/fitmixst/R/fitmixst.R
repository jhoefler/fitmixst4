# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Fitted mixed model
#' @S3method fitmixst numeric
#' @param y a multidimensional input vector.
#' @param ... the rest
#' @return a object of the class fitmixst.
#' @keywords fit mixed skew t.
#' @export
#' @examples
#' 
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)
#' out <- fitmixst(y,g=2,method="kmeans")

fitmixst <- function(y, ...) UseMethod("fitmixst")

