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

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)

test_cluster_analysis <- function() {
    checkException(cluster_analysis(sel.exp=test_matrix, cluster_type=1,
        distance="euclidean", linkage_type="ward.D2",
        gene_distance="correlation", num_clusters=3, data_name="word",
        probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
        cluster_num_selection="Fixed_Clust_Num"))
}
