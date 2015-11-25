######## Functions to Normalize Data and Write Matrices to Text Files ###############

#' Function to write a data matrix to a text file.
#' @param	 tmpMatrix The object matrix containing data.
#' @param tmpFileName The string name of the text file to write the matrix to.
#' @param blnRowNames Logical value indicating if row names of the matrix should be written
#' along with the matrix.
#' @param blnColNames Logical value indicating if the column names of the matrix should be
#' written with the matrix.
#' @return Text file containing the data matrix.
#' @author Peiyong Guan
#' @seealso \code{\link[utils]{write.table}}
#' @examples
#'
#' #Load sample dataset
#' data(iris)
#' # View sample matrix
#' iris[1:4,1:4]
#' # Write sample matrix to text file
#' WriteMatrixToFile(tmpMatrix = iris, tmpFileName = "iris.sample.matrix.txt",
#' blnRowNames = TRUE, blnColNames = TRUE)
#'
#' @export
WriteMatrixToFile <- function(tmpMatrix, tmpFileName, blnRowNames, blnColNames)
{
  output <- file(tmpFileName, "at")
  utils::write.table(tmpMatrix, output, sep = "\t", quote = FALSE,
              row.names = blnRowNames, col.names = blnColNames)
  close(output)
}

