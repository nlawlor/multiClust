number_clusters <- function(data.exp, Fixed=3, gap_statistic=NULL){
  # Conditional to cause error if user chooses both Fixed and
  #gap_statistic method
  if (is.null(Fixed) == FALSE) {
    if (is.null(gap_statistic) == FALSE) {
      stop("Please choose one of the two options
           (Fixed or gap_statistic) only")
    }
  }

  if (is.null(Fixed) == FALSE) {
    if (is.numeric(Fixed) == TRUE) {
      cluster_num <- Fixed
      print(paste("The fixed cluster number is: ", cluster_num, sep=""))
      return(cluster_num)
    } else {
      stop("Please input an integer for Fixed")
    }
  }
  if (gap_statistic == TRUE) {
    #Makes Expression Data Numeric
    colname <- colnames(data.exp)
    data.exp <- t(apply(data.exp, 1, as.numeric))
    colnames(data.exp) <- colname

    # Transpose the matrix to do the gap_stat calculation by columns
    # clusGap function calculates by rows
    theGap <- cluster::clusGap(t(data.exp),
       FUNcluster=stats::kmeans, K.max=8)
    # This returns a table with K.max rows and 4 other columns
    # The maximum "gap" value corresponds to the optimal cluster number
    gapDF <- as.data.frame(theGap$Tab)
    # Make a new data frame of cluster number and gap statistic
    clust_total <- c(1:dim(gapDF)[1])
    gap_and_clust <- data.frame(clust_total, gapDF[, "gap"])
    # Determine which is maximum gap statistic
    max_gap <- max(gap_and_clust[,2])
    # Obtain cluster number with maximum gap statistic
    cluster_num2 <- which(gap_and_clust[,2] == max_gap)
    print(paste("The gap statistic cluster number is: ",
       cluster_num2, sep=""))
    return(cluster_num2)
  }
}

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)

test_number_clusters <- function() {
    checkException(number_clusters(data.exp=test_matrix, Fixed="string",
        gap_statistic=NULL), msg="Unable to use a string")
}

