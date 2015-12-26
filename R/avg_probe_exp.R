## VI. Average Expression of Each Gene Probe in Each Sample Cluster ##
#' Function to produce a matrix containing the average expression of each
#' gene probe within each sample cluster.
#' @param sel.exp Object containing the numeric selected gene
#' expression matrix. This object is an output of the probe_ranking function.
#' @param samp_cluster Object vector containing the samples and the cluster
#' number they belong to. This is an output of the cluster_analysis function.
#' @param data_name String indicating the cancer type and name of the
#' dataset being analyzed. This name will be used to label the sample
#' dendrograms and heatmap files.
#' @param cluster_type String indicating the type of clustering method
#' used in the cluster_analysis function. "Kmeans" or "HClust"
#' are the two options.
#' @param distance String describing the distance metric uses for HClust in
#' the cluster_analysis function. Options include one of "euclidean",
#' "maximum", manhattan", "canberra", "binary", or "minkowski".
#' @param linkage_type String describing the linkage metric used in the
#' cluster_analysis function. Options include "ward.D2", "average",
#' "complete", "median", "centroid", "single", and "mcquitty".
#' @param probe_rank String indicating the feature selection method used
#' in the probe_ranking function. Options include "CV_Rank", "CV_Guided",
#' "SD_Rank", and "Poly".
#' @param probe_num_selection String indicating the way in which probes
#' were selected in the number_probes function. Options include
#' "Fixed_Probe_Num", "Percent_Probe_Num", and "Adaptive_Probe_Num".
#' @param cluster_num_selection String indicating how the number of
#' clusters were determined in the number_clusters function. Options include
#' "Fixed_Clust_Num" and "Gap_Statistic".
#' @author Nathan Lawlor, Alec Fabbri
#' @seealso \code{\link{number_clusters}}, \code{\link{number_probes}},
#' \code{\link{probe_ranking}}, \code{\link{cluster_analysis}}
#' @return Returns an object matrix with the average mean expression for
#' each probe in each sample cluster. Also outputs the object to a text file.
#' @examples
#' # Produce matrix of average expression of each probe in each cluster
#' # Load in a data file
#' data_file <- system.file("extdata", "GSE2034.normalized.expression.txt",
#'     package="multiClust")
#' data <- input_file(input=data_file)
#' # Choose 300 genes to select for
#' gene_num <- number_probes(input=data_file, data.exp=data, Fixed=300,
#'    Percent=NULL, Adaptive=NULL)
#' # Choose the "CV_Rank" Method for gene ranking
#' sel.data <- probe_ranking(input=data_file, probe_number=300,
#'     probe_num_selection="Fixed_Probe_Num", data.exp=data, method="CV_Rank")
#' # Choose a fixed cluster number of 3
#' clust_num <- number_clusters(data.exp=data, Fixed=3, gap_statistic=NULL)
#' # Call function for Kmeans parameters
#' kmeans_analysis <- cluster_analysis(sel.exp=sel.data, cluster_type="Kmeans",
#'     distance=NULL, linkage_type=NULL, gene_distance=NULL,
#'     num_clusters=3, data_name="GSE2034 Breast",
#'     probe_rank="CV_Rank", probe_num_selection="Fixed_Probe_Num",
#'     cluster_num_selection="Fixed_Clust_Num")
#' # Call function for average matrix expression calculation
#' avg_matrix <- avg_probe_exp(sel.exp=sel.data, samp_cluster=kmeans_analysis,
#'     data_name="GSE2034 Breast", cluster_type="Kmeans", distance=NULL,
#'     linkage_type=NULL, probe_rank="CV_Rank",
#'     probe_num_selection="Fixed", cluster_num_selection="Fixed_Clust_Num")
#' @export

avg_probe_exp <- function(sel.exp, samp_cluster, data_name,
    cluster_type="HClust", distance="euclidean", linkage_type="ward.D2",
    probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num") {

    # Conditionals to make sure inputs are strings
    if (is.character(cluster_type) == FALSE) {
        stop("Please input string for cluster_type")
    }

    if (is.null(distance) == FALSE) {
        if (is.character(distance) == FALSE) {
            stop("Please input string for distance")
        }
    }

    if (is.null(linkage_type) == FALSE) {
        if (is.character(linkage_type) == FALSE) {
            stop("Please input string for linkage_type")
        }
    }

    if (is.character(data_name) == FALSE) {
        stop("Please input string for data_name")
    }

    if (is.character(probe_rank) == FALSE) {
        stop("Please input string for probe_rank")
    }

    if (is.character(probe_num_selection) == FALSE) {
        stop("Please input string for probe_num_selection")
    }

    if (is.character(cluster_num_selection) == FALSE) {
        stop("Please input string for cluster_num_selection")
    }

    # Function to write matrix to text file
    WriteMatrixToFile <- function(tmpMatrix, tmpFileName,
        blnRowNames, blnColNames) {
        output <- file(tmpFileName, "at")
        utils::write.table(tmpMatrix, output, sep="\t", quote=FALSE,
            row.names=blnRowNames, col.names=blnColNames)
        close(output)
    }

    # Make gene matrix numeric
    colname <- colnames(sel.exp)
    sel.exp <- t(apply(sel.exp, 1, as.numeric))
    colnames(sel.exp) <- colname

    # Maximum number of clusters
    numofclust <- max(samp_cluster)
    # Produce a matrix with average gene expression in each cluster
    meansofclust <- NULL
    names <- NULL

    for (i in 1:numofclust) {
        # names contains the column names specifying the cluster number
        names <- c(names, paste('Cluster', i, sep=' '))
        # savednames contains samples that belong to a given cluster number
        savednames <- which(samp_cluster == i)
        vec <- NULL
        for (i in 1:length(savednames)) {
            # match all the sample names in the selected exp file
            # with sample names from savednames
            num <- which(colnames(sel.exp) == names(savednames)[i])
            # Vec will hold all indices of the sample names
            vec <- c(vec, num)
        }
        snipexp = sel.exp[,vec]
        # If statement to catch if only one sample in the cluster
        if (is.null(dim(snipexp)[1]) == TRUE) {
            row.mean <- snipexp
        } else {
            row.mean <- apply(snipexp, 1, mean)
        }
        meansofclust <- rbind(meansofclust, row.mean)
    }

    meansofclust <- t(meansofclust)
    colnames(meansofclust) <- names
    # Write matrix to file
    WriteMatrixToFile(meansofclust, paste(data_name, cluster_type, distance,
        linkage_type, probe_rank, probe_num_selection,
        cluster_num_selection, "Avg.Exp.txt", sep="."), TRUE, TRUE)
    print(paste("Your matrix containing the average gene probe expression",
        "in each cluster is finished"))
    return(meansofclust)
}
