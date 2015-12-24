############# V. Clustering Option ##############
#' Function to perform Kmeans or Hierarchical clustering analysis of the
#' selected gene probe expression data.
#' @param sel.exp Object containing the numeric selected gene expression
#' matrix. This object is an output of the probe_ranking function.
#' @param cluster_type String indicating the type of clustering method to use.
#' "Kmeans" or "HClust" are the two options. The default is set to "HClust".
#' @param distance String describing the distance metric to use for
#' the dist function during hierarchical clustering. dist uses a default
#' distance metric of Euclidean distance. Options include one of "euclidean",
#' "maximum", manhattan", "canberra", "binary", or "minkowski". Kmeans
#' clustering does not use a distance metric. The default value is set
#' to "euclidean".
#' @param linkage_type String describing the linkage metric to be
#' used for HClust. The default is set to "ward.D2", however other options
#' include "average", "complete", "median", "centroid",
#' "single", and "mcquitty".
#' @param gene_distance String describing the distance measure to be used for
#' the Dist function when performing hierarchical clustering of genes.
#' Options include one of "euclidean", "maximum", "manhattan", "canberra",
#' "binary", "pearson", "abspearson", "correlation", "abscorrelation",
#' "spearman" or "kendall". The default of gene_distance is set
#' to "correlation". The deafult value is set to "correlation".
#' The argument can be set to NULL when Kmeans clustering is used.
#' @param num_clusters Positive integer to specify the number of clusters
#' samples will be divided into. This number is determined by the
#' number_clusters function.
#' @param data_name String indicating the cancer type and name of the
#' dataset being analyzed. This name will be used to label the sample
#' dendrograms and heatmap files.
#' @param probe_rank String indicating the feature selection method used
#' in the probe_ranking function. Options include "CV_Rank", "CV_Guided",
#' "SD_Rank", and "Poly".
#' @param probe_num_selection String indicating the way in which probes
#' were selected in the number_probes function. Options include
#' "Fixed_Probe_Num", "Percent_Probe_Num", and "Adaptive_Probe_Num".
#' @param cluster_num_selection String indicating how the number of clusters
#' were determined in the number_clusters function. Options include
#' "Fixed_Clust_Num" and "Gap_Statistic".
#' @return Returns a vector containing the sample information and respective
#' cluster number. In addition, this function outpus sample cluster
#' dendrogams, average expression for each probe in each cluster, and
#' heatmap images and Java TreeView files for HClust dendrograms.
#' @seealso \code{\link{probe_ranking}}, \code{\link{number_clusters}},
#' \code{\link{number_probes}}, \code{\link[stats]{hclust}},
#' \code{\link[stats]{kmeans}}, \code{\link[amap]{Dist}},
#' \code{\link[stats]{dist}}
#' @author Nathan Lawlor, Alec Fabbri
#' @examples
#'
#' # Example 1: HClust Analysis
#' # Load in a data file
#' data_file <- system.file("extdata", "GSE2034.normalized.expression.txt",
#'     package="multiClust")
#' data <- input_file(input=data_file)
#' # Choose 300 genes to select for
#' gene_num <- number_probes(input=data_file, data.exp=data, Fixed=300,
#'     Percent=NULL, Adaptive=NULL)
#' # Choose the "CV_Rank" Method for gene ranking
#' sel.data <- probe_ranking(input=data_file, probe_number=300,
#'     probe_num_selection="Fixed_Probe_Num", data.exp=data, method="CV_Rank")
#' # Choose a fixed cluster number of 3
#' clust_num <- number_clusters(data.exp=data, Fixed=3, gap_statistic=NULL)
#'
#' # Call function using HClust parameters
#' hclust_analysis <- cluster_analysis(sel.exp=sel.data, cluster_type="HClust",
#'     distance="euclidean", linkage_type="ward.D2",
#'     gene_distance="correlation", num_clusters=3,
#'     data_name="GSE2034 Breast", probe_rank="CV_Rank",
#'     probe_num_selection="Fixed_Probe_Num",
#'     cluster_num_selection="Fixed_Clust_Num")
#'
#' # Example 2: Kmeans Analysis
#' # Call function for Kmeans parameters
#' kmeans_analysis <- cluster_analysis(sel.exp=sel.data, cluster_type="Kmeans",
#'     distance=NULL, linkage_type=NULL, gene_distance=NULL,
#'     num_clusters=3, data_name="GSE2034 Breast",
#'     probe_rank="CV_Rank", probe_num_selection="Fixed_Probe_Num",
#'     cluster_num_selection="Fixed_Clust_Num")
#' @export
cluster_analysis <- function(sel.exp, cluster_type="HClust",
    distance="euclidean", linkage_type="ward.D2",
    gene_distance="correlation", num_clusters, data_name,
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

    if (is.null(gene_distance) == FALSE) {
        if (is.character(gene_distance) == FALSE) {
            stop("Please input string for gene_distance")
        }
    }

    if (is.numeric(num_clusters) == FALSE) {
        stop("Please input integer for num_clusters")
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

    # Function to normalize the gene expression matrix
    nor.min.max <- function(x){
        x.min <- min(x)
        x.max <- max(x)
        x <- (x-x.min) / (x.max - x.min)
        return (x)
    }

    # Ensure the selected gene expression matrix is numeric
    colname <- colnames(sel.exp)
    sel.exp <- t(apply(sel.exp, 1, as.numeric))
    colnames(sel.exp) <- colname

    # Normalization of the selected gene expression matrix
    exp.sel <- t(apply(sel.exp, 1, nor.min.max))
    # Code for HClust analysis of data
    if (cluster_type == "HClust") {
        # Hclust analysis, divide samples into cluster number specified
        hc <- stats::hclust(stats::dist(x=t(exp.sel), method=distance),
            method=linkage_type)
        dc <- stats::as.dendrogram(hc)
        op <- graphics::par(mar=c(1,2,2,1) + 0.1)
        dc <- dendextend::set(dc, "branches_k_color", value=1:num_clusters,
            k=num_clusters)

        # Make a vector filled with blank spaces to replace sample
            # names with blank space
        vec <- NULL
        len <- length(labels(dc))
        for (i in 1:len) {
            i <- ""
            vec <- c(vec, i)
        }

        # Make another dendrogram with no labels
        dc1 <- dc
        # Remove labels from this dendrogram
        dc1 <- dendextend::set(dc1, "labels", vec)

        # Output the sample dendrogram to a pdf
        grDevices::pdf(paste(data_name, distance, linkage_type, probe_rank,
            probe_num_selection , cluster_num_selection, "pdf", sep="."))
        graphics::plot(dc1, main=paste(data_name, "Sample Dendrogram",
            sep=" "), ylab="Height")
        graphics::par(op)
        dendextend::rect.dendrogram(dc1, k=num_clusters, border='orange')
        grDevices::dev.off()
        print("Your HClust Sample Dendrogram has been outputted")

        # Output heatmap files to view in Java TreeView
        # Take the mean difference of the selected gene expression matrix
        exp.sel.mean <- apply(sel.exp, MARGIN=1, mean)
        new.exp.sel <- sel.exp - exp.sel.mean

        # Obtain the sample tree file
        hcSample <- stats::hclust(stats::dist(x=t(new.exp.sel),
            method=distance), method=linkage_type)
        ctc::r2atr(hcSample, file=paste(data_name, probe_rank, linkage_type,
            probe_num_selection, cluster_num_selection, "atr", sep="."),
            distance=hcSample$method)

        # Output the gene tree file
        hcGene <- stats::hclust(amap::Dist(x=new.exp.sel,
            method=gene_distance), method=linkage_type)
        ctc::r2gtr(hcGene, file=paste(data_name, probe_rank, linkage_type,
            probe_num_selection, cluster_num_selection, "gtr", sep="."),
            distance=hcGene$method)

        # Output the CDT tree file
        ctc::r2cdt(hcGene, hcSample, new.exp.sel, file=paste(data_name,
            probe_rank, linkage_type, probe_num_selection,
            cluster_num_selection, "cdt", sep="."))
        print(paste("Your atr, gtr, and cdt files have been outputted",
                    "for viewing in Java TreeView"))

        # Divide the samples into respective clusters
        group <- dendextend::cutree(dc, num_clusters,
            order_clusters_as_data=FALSE)
        group <- group[order(names(group))]
        utils::write.csv(group, paste(data_name, cluster_type, distance,
            linkage_type, probe_rank, probe_num_selection,
            cluster_num_selection, "Samples.Clusters.csv"))
        print(paste("A CSV file has been produced containing your sample",
            "and cluster assignment information"))
        # Return sample and cluster information
        return(group)
    }

    # Code for Kmeans clustering analysis
    if (cluster_type == "Kmeans"){
        # Kmeans clustering and plot output
        set.seed(1293075)
        exp.sel.k <- stats::kmeans(t(exp.sel), centers=num_clusters, nstart=25)
        kmeans_name <- paste(data_name, probe_rank, probe_num_selection,
            cluster_num_selection,"Kmeans.Plot.pdf", sep = ".")

        # Obtain vector showing which  cluster each sample belongs to
        samp_clusters <- exp.sel.k$cluster
        utils::write.csv(samp_clusters, paste(data_name, cluster_type,
            probe_rank, probe_num_selection, cluster_num_selection,
            "Samples.Clusters.csv"))
        print(paste("A CSV file has been produced containing your sample",
            "and cluster assignment information"))
        return(samp_clusters)
    }
}
