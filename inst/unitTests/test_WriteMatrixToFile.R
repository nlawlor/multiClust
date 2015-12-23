WriteMatrixToFile <- function(tmpMatrix, tmpFileName, blnRowNames,
    blnColNames) {
    output <- file(tmpFileName, "at")
    write.table(tmpMatrix, output, sep="\t", quote=FALSE,
        row.names=blnRowNames, col.names=blnColNames)
    close(output)
}

test_WriteMatrixToFile <- function() {
    checkException(WriteMatrixToFile(tmpMatrix= "string",
        tmpFileName= 1, blnRowNames=TRUE,
        blnColNames=TRUE), msg="Unable to use a number for file name")
}
