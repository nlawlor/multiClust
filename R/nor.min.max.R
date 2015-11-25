#' Function to normalize data to bring values into alignment. This function
#' uses feature scaling to normalize values in a dataset between 0 and 1.
#' @param x An integer object of numeric value
#' @return Returns a numeric value normalized between 0 and 1.
#' @author Peiyong Guan
#' @examples
#'
#' # Load sample dataset
#' data(iris)
#' # View sample matrix
#' iris[1:5,1:5]
#' # Coerce sample matrix to numeric values
#' iris <- t(apply(iris[,1:4], 1, as.numeric))
#' #Normalize values in the matrix using the function
#' data.min.max <- t(apply(iris, 1, nor.min.max))
#' @export
nor.min.max <- function(x) {
  x.min <- min(x)
  x.max <- max(x)
  x <- (x - x.min) / (x.max - x.min)
  return (x)
}
