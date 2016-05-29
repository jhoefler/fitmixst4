#' Write matrix to console to copy paste in latex
#' 
#' @param mat matrix 
#' @param digits round to those digits. Defaults to 2.
#' @return NULL
#' @description Just prints in console
#' @author Andi Boeck
#' @export
#' @examples
#'  mat <- matrix(rnorm(100), 10,10)
#'  toLatexMat(mat)

toLatexMat <- function(mat, digits=2)
{
  mat <- round(mat,digits)
  tmat <- t(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  cat("\\left[", "\n")
  cat(paste("\\begin{array}", "{", 
            paste(rep("c", nc), sep= "", collapse = ""), "}", collapse= "", sep=""), "\n")
  for (i in 1:nr)
  {
    cat(paste(paste(mat[i,], collapse= " & "), "\\\\"), "\n")
  }
  cat("\\end{array}", "\n")
  cat("\\right]", "\n")
}

