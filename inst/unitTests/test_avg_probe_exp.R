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

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)
samples <- c("a", "b", "c,", "d")

test_avg_probe_exp <- function() {
    checkException(avg_probe_exp(sel.exp=test_matrix, samp_cluster=samples,
    cluster_type=1, distance="euclidean", linkage_type="ward.D2",
    data_name="word", probe_rank="SD_Rank",
    probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num"))
}

