######## I. Input Gene Expression Matrix ####################
#' Function to read-in the gene expression file and assign gene probe names as the rownames.
#' @param input String indicating the name of the text file containing the gene
#' expression matrix to be read in. This matrix file should have the gene probes
#' in the first column of the matrix. The gene probes will be assigned as
#' the rownames of the matrix.
#' @return Returns an object containing the gene expression matrix with the
#' gene probe names as the rownames.
#' @seealso \code{\link[utils]{read.table}}
#' @author Nathan Lawlor
#' @note This function works best when using gene expression datasets from Gene Expression
#' Omnibus.
#' @examples
#' # Load in a test file
#' data_file <- system.file("extdata", "GSE2034.normalized.expression.txt", package = "multiClust")
#' data <- input_file(input = data_file)
#' # View matrix with gene probes assigned as rownames
#' data[1:4,1:4]
#' @export
input_file <- function(input) {
  data.exp <- utils::read.delim2(input, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)
  rownames(data.exp) <- data.exp[,1]
  data.exp[,1] <- NULL
  print("The gene expression matrix has been loaded")
  return(data.exp)
}
