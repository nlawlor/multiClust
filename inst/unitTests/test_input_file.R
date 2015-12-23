input_file <- function(input) {
    data.exp <- utils::read.delim2(input, header=TRUE,
        stringsAsFactors=FALSE, row.names=NULL)
    rownames(data.exp) <- data.exp[, 1]
    data.exp[, 1] <- NULL
    print("The gene expression matrix has been loaded")
    return(data.exp)
}

test_input_file <- function() {
  checkException(input_file(input=1), msg="Unable to use integer as argument")
}
