#' refit multivariate linear model
#' 
#' @param mod suitable model
#' @return An object of class \code{"mle2"} \link{mle2-class}
#' @seealso mle2
#' @seealso mle2-class
#' @import bbmle
#' @examples 
#' mod <- lm(cbind(Fertility, Agriculture, Infant.Mortality)~ 
#'  Examination + Education  + Catholic, data=swiss)
#' start <- get.theta(mod)
#' fit <- refit.mlm(mod, start=start)
#' pr <- profile(fit, which=13) # sigma2 of first response column
#' plot(pr)

refit.mlm <-
function(mod, start=NULL, ...)
{
  ####### define function with seperate arguments
  arglist <- get.theta(mod)
  args <- paste(names(arglist), collapse=",")
  string3 <- paste("mlmloglik(c(", args, "))", sep="")
  exp1 <- parse(text=string3)
  multipar <- function(x) a+b
  formals(multipar) <- arglist
  body(multipar) <- as.call(c(as.name("{"), exp1))  
  ######################################
  
  
  fit <- mle2(multipar, start=start,...)
  
  return(fit)
  
}
