surv_analysis <- function(samp_cluster, clinical, survival_type="RFS",
    data_name, cluster_type="HClust", distance="euclidean",
    linkage_type="ward.D2", probe_rank="SD_Rank",
    probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num") {

  # Conditionals to make sure inputs are strings
  if (is.character(clinical) == FALSE) {
    stop("Please input string for clinical")
  }

  if (is.character(survival_type) == FALSE) {
    stop("Please input string for survival_type")
  }

  if (is.character(data_name) == FALSE) {
    stop("Please input string for data_name")
  }

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

  if (is.character(probe_rank) == FALSE) {
    stop("Please input string for probe_rank")
  }

  if (is.character(probe_num_selection) == FALSE) {
    stop("Please input string for probe_num_selection")
  }

  if (is.character(cluster_num_selection) == FALSE) {
    stop("Please input string for cluster_num_selection")
  }

  # Read in survival data text file
  sample.anns <- utils::read.delim2(clinical, header=TRUE,
                                    stringsAsFactors=FALSE)
  time <- sample.anns[, 1]
  event <- sample.anns[, 2]
  time <- as.numeric(time)
  event <- as.numeric(event)
  event_type <- survival_type
  time_type <- 'Months'

  # Determine max number of clusters in samp_cluster
  Number_Clusters <- max(samp_cluster)

  # Produce Kaplan Meier survival plots
  lev <- '1'
  for(i in 2:Number_Clusters) {
    lev <- c(lev, paste('',i, sep=''))
  }

  groupfac <- factor(samp_cluster, levels=lev)
  cox <- survival::coxph(survival::Surv(time, event==1) ~ groupfac)
  csumm <- summary(cox)
  cox.p.value <- c(csumm$logtest[3])
  p <- survival::survfit(survival::Surv(time,
    event==1)~survival::strata(samp_cluster))

  # Produce Kaplan-Meier Survival Plot
  grDevices::pdf(paste(data_name, cluster_type, linkage_type, probe_rank,
    probe_num_selection, cluster_num_selection,
    survival_type, "pdf", sep='.'))
  graphics::plot(p, xlab= time_type, ylab='Survival Probability',
    conf.int=FALSE, xlim=c(0,100), col=1:Number_Clusters,
    main=paste(data_name, cluster_type, probe_rank,
    survival_type, sep =' '))
  graphics::legend('bottomleft',
    legend=levels(survival::strata(samp_cluster)),lty=1,
    col=1:Number_Clusters,text.col=1:Number_Clusters, title='clusters')
  graphics::legend('topright', paste('Pvalue =',
    signif(csumm$logtest[3],digits=2), sep=''))
  grDevices::dev.off()
  print("Your Kaplan Meier Survival Plot has been finished")
  return(cox.p.value)
}

# Make a test samp_cluster vector
samples <- c("a", "b", "c,", "d")

test_surv_analysis <- function() {
    checkException(surv_analysis(samp_cluster=samples, clinical = 1,
    surv_type = "RFS", data_name = "test", cluster_type="HClust",
    distance="euclidean", linkage_type="ward.D2",
    data_name="word", probe_rank="SD_Rank",
    probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num"))
}
