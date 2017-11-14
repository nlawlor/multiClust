nor.min.max <- function(x) {
    if (is.numeric(x) == FALSE) {
        stop("Please input numeric for x")
    }
    x.min <- min(x)
    x.max <- max(x)
    x <- (x - x.min) / (x.max - x.min)
    return (x)
}

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)

# Make a normalized matrix
normalized_matrix <- nor.min.max(test_matrix)

test_nor.min.max <- function() {
    checkEqualsNumeric(nor.min.max(test_matrix), normalized_matrix)
    checkException(nor.min.max("a"), msg="Unable to use a string")
}
