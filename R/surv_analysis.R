#############  VII. Survival Analysis  ###################
#' Function to produce Kaplan-Meier Survival Plots of selected gene expression data.
#' @param samp_cluster Object vector containing the samples and the cluster number they belong to.
#' This object is an output of the cluster_analysis function.
#' @param clinical String indicating the name of the text file containing patient clinical information.
#' The file should be a data frame consisting of two columns. The first column contains the
#' patient survival time information in months. The second column indicates occurrence of
#' a censorship (0) or an event (1).
#' @param survival_type String specifying the type of survival event being analyzed. Examples include
#' "Disease-free survival (DFS)", "Overall Survival (OS)", "Relapse-free survival (RFS)", etc.
#' @param data_name String indicating the name to be used to label the plot.
#' @param cluster_type String indicating the type of clustering method used in the
#' cluster_analysis function. "Kmeans" or "HClust" are the two options.
#' @param distance String describing the distance metric uses for HClust in
#' the cluster_analysis function. Options include one of "euclidean", "maximum", manhattan",
#' "canberra", "binary", or "minkowski".
#' @param linkage_type String describing the linkage metric used in the
#' cluster_analysis function. Options include "ward.D2", "average", "complete",
#' "median", "centroid", "single", and "mcquitty".
#' @param probe_rank String indicating the feature selection method used in the
#' probe_ranking function. Options include "CV_Rank", "CV_Guided", "SD_Rank",
#' and "Poly".
#' @param probe_num_selection String indicating the way in which probes were selected
#' in the number_probes function. Options include "Fixed_Probe_Num",
#' "Percent_Probe_Num", and "Adaptive_Probee_Num".
#' @param cluster_num_selection String indicating how the number of clusters were
#' determined in the number_clusters function.Options include "Fixed_Clust_Num"
#'  and "Gap_Statistic".
#' @return Produces a pdf image of a Kaplan-Meier Survival Plot with Cox Survival P Value.
#' Also returns an object containing the cox survival P value.
#' @seealso \code{\link{number_clusters}}, \code{\link{number_probes}},
#' \code{\link{probe_ranking}}, \code{\link{cluster_analysis}},
#' \code{\link[survival]{coxph}}
#' @author Alec Fabbri, Nathan Lawlor
#' @examples
#' # Load in a data file
#' data_file <- system.file("extdata", "GSE2034.normalized.expression.txt", package = "multiClust")
#' data <- input_file(input = data_file)
#' # Choose 300 genes to select for
#' gene_num <- number_probes(input = data, data.exp = data, Fixed = 300,
#'    Percent = NULL, Adaptive = NULL)
#' # Choose the "CV_Rank" Method for gene ranking
#' sel.data <- probe_ranking(input = data_file, probe_number = 300,
#'    probe_num_selection = "Fixed_Probe_Num", data.exp = data, method = "CV_Rank")
#' # Choose a fixed cluster number of 3
#' clust_num <- number_clusters(data.exp = data, Fixed = 3, gap_statistic = NULL)
#' # Call function for Kmeans parameters
#' kmeans_analysis <- cluster_analysis(sel.exp = sel.data, cluster_type = "Kmeans",
#'    distance = NULL, linkage_type = NULL, gene_distance = NULL,
#'    num_clusters = 3, data_name = "GSE2034 Breast",
#'    probe_rank = "CV_Rank", probe_num_selection = "Fixed_Probe_Num",
#'    cluster_num_selection = "Fixed_Clust_Num")
#' # Load the clinical outcome file
#' clin_file <- system.file("extdata", "GSE2034-RFS-clinical-outcome.txt", package = "multiClust")
#' # Example of Calling surv_analysis function
#' surv <- surv_analysis(samp_cluster = kmeans_analysis, clinical = clin_file,
#'    survival_type = "RFS", data_name = "GSE2034 Breast", cluster_type = "Kmeans",
#'    distance = NULL, linkage_type = NULL, probe_rank = "CV_Rank",
#'    probe_num_selection = "Fixed_Probe_Num",
#'    cluster_num_selection = "Fixed_Cluster_Num")
#'
#' @export
surv_analysis <- function(samp_cluster, clinical, survival_type = "RFS", data_name,
                          cluster_type = "HClust", distance = "euclidean",
                          linkage_type = "ward.D2", probe_rank = "SD_Rank",
                          probe_num_selection = "Fixed_Probe_Num",
                          cluster_num_selection = "Fixed_Clust_Num"){
  # Read in survival data text file
  sample.anns <- utils::read.delim2(clinical, header = TRUE, stringsAsFactors = FALSE)
  time <- sample.anns[,1]
  event <- sample.anns[,2]
  time <- as.numeric(time)
  event <- as.numeric(event)
  event_type <- survival_type
  time_type <- 'Months'

  # Determine max number of clusters in samp_cluster
  Number_Clusters <- max(samp_cluster)

  # Produce Kaplan Meier survival plots
  lev <- '1'
  for(i in 2:Number_Clusters) {
    lev <- c(lev, paste('',i, sep = ''))
  }

  groupfac <- factor(samp_cluster, levels = lev)
  cox <- survival::coxph(survival::Surv(time, event==1) ~ groupfac)
  csumm <- summary(cox)
  cox.p.value <- c(csumm$logtest[3])
  p <- survival::survfit(survival::Surv(time,event==1)~survival::strata(samp_cluster))

  # Produce Kaplan-Meier Survival Plot
  grDevices::pdf(paste(data_name, cluster_type, linkage_type, probe_rank, probe_num_selection, cluster_num_selection, survival_type, "pdf", sep = '.'))
  graphics::plot(p,xlab= time_type ,ylab='Survival Probability',conf.int=FALSE, xlim = c(0,100), col=1:Number_Clusters, main = paste(data_name, cluster_type, probe_rank, survival_type, sep =' '))
  graphics::legend('bottomleft',legend=levels(survival::strata(samp_cluster)),lty=1,col=1:Number_Clusters,text.col=1:Number_Clusters, title='clusters')
  graphics::legend('topright', paste('Pvalue =', signif(csumm$logtest[3],digits = 2), sep =''))
  grDevices::dev.off()

  print("Your Kaplan Meier Survival Plot has been finished")
  return(cox.p.value)
}
