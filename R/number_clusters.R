########## IV. Number of Clusters  ######################
#' Function to determine the number of clusters to be used to cluster gene probes and samples.
#' @param data.exp The numeric original gene expression matrix to be used for clustering of
#' genes and samples. This object is an output of the input_file function.
#' @param Fixed A positive integer used to represent the number of clusters the samples
#' and probes will be divided into.
#' @param gap_statistic A logical indicating whether to use the gap_statistic to
#' determine the optimal number of clusters to divide samples into.
#' @note The user should only choose either the fixed or gap_statistic option, not both.
#' When using the gap_statistic option, change the argument to TRUE and "Fixed" to NULL.
#' @return An object with the determined number of clusters to use.
#' @author Nathan Lawlor, Alec Fabbri
#' @seealso \code{\link[cluster]{clusGap}}, \code{\link{probe_ranking}}
#' @examples
#'
#' #Example 1: Using a fixed cluster number
#' # Load in a test file
#' data_file <- system.file("extdata", "GSE2034.normalized.expression.txt", package = "multiClust")
#' data <- input_file(data_file)
#' clust_num <- number_clusters(data.exp = data, Fixed = 3, gap_statistic = NULL)
#'
#' \dontrun{
#' # Example 2: Using the gap_statistic to determine the optimal cluster number
#' # Computation time is somewhat long
#' clust_num <- number_clusters(data.exp = data, Fixed = NULL, gap_statistic = TRUE)
#'  }
#' @export
number_clusters <- function(data.exp, Fixed = 3, gap_statistic = NULL){
  # Conditional to cause error if user chooses both Fixed and gap_statistic method
  if (is.null(Fixed) == FALSE) {
    if (is.null(gap_statistic) == FALSE) {
      stop("Please choose one of the two options (Fixed or gap_statistic) only")
    }
  }

  if (is.null(Fixed) == FALSE) {
    cluster_num = Fixed
    print(paste("The fixed cluster number is: ", cluster_num, sep = ""))
    return(cluster_num)
  }

  if (gap_statistic == TRUE) {

    #Makes Expression Data Numeric
    colname = colnames(data.exp)
    data.exp = t(apply(data.exp, 1,as.numeric))
    colnames(data.exp) = colname

    # Transpose the matrix to do the gap_stat calculation by columns (samples)
    # clusGap function calculates by rows
    theGap <- cluster::clusGap(t(data.exp), FUNcluster = stats::kmeans, K.max = 8)
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
    print(paste("The gap statistic cluster number is: ", cluster_num2, sep = ""))
    return(cluster_num2)
  }

}
