number_probes <- function(input, data.exp, Fixed=1000, Percent=NULL,
    Poly=NULL, Adaptive=NULL) {
  # Conditional to cause error if user chose multiple methods at once
  if (is.null(Fixed) == FALSE) {
    if (is.null(Percent) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
    }

  if (is.null(Fixed) == FALSE) {
    if (is.null(Adaptive) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
    }

  if (is.null(Percent) == FALSE) {
    if (is.null(Adaptive) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
    }

  if (is.null(Fixed) == FALSE) {
    if (is.null(Poly) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
    }

  if (is.null(Percent) == FALSE) {
    if (is.null(Poly) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
    }

  if (is.null(Adaptive) == FALSE) {
    if (is.null(Poly) == FALSE) {
      stop("Please choose one of the four options
           (Fixed, Percent, Poly, Adaptive) only")
    }
  }

  # For Fixed method
  if (is.null(Fixed) == FALSE) {
    # Conditional to make sure Fixed is numeric
    if (is.numeric(Fixed) == FALSE) {
      stop("Please input numeric for Fixed")
    }
    gene_num=Fixed
    # Conditional to produce error if fixed gene number is
    # greater than total genes in matrix
    if (gene_num > dim(data.exp)[1]) {
      stop("Gene number chosen is too large")
    }
    print(paste("The fixed gene probe number is: ", gene_num, sep=""))
    return(gene_num)
  }

  # For percent method
  if (is.null(Percent) == FALSE) {
    # Conditional to produce error if Percent chosen is not greater
    # than 0 or less/equal to 1
    if (Percent < 0) {
      stop("Please choose a positive number between 0 and 100")
    }
    if (Percent > 100) {
      stop ("Please choose a positive number between 0 and 100")
    }
    gene_num2 <- round((Percent / 100) * dim(data.exp)[1])
    print(paste("The percent gene probe number is: ", gene_num2, sep=""))
    return(gene_num2)
  }

  # For Poly method
  if (is.null(Poly) == FALSE) {
    # Make data matrix numeric
    colname <- colnames(data.exp)
    data.exp <- t(apply(data.exp, 1, as.numeric))
    colnames(data.exp) <- colname

    # Get mean and standard deviation of each gene
    row.mean <- apply(data.exp, 1, mean)
    row.sd <- apply(data.exp, 1, stats::sd)
    transinput <- data.exp

    # Fit a second degree polynomial to data
    fit <- stats::lm(row.sd ~ poly(row.mean, 2, raw=TRUE))
    res <- fit$residuals
    # Determine which genes are above the fitted polynomial
    numgenesel <- length(which(res > 0))
    selgenes <- names(which(res > 0))

    # Obtain the selected gene matrix for genes above the polynomial
    vec <- NULL
    for (i in 1:length(selgenes)){
      num <- which(rownames(transinput) == selgenes[i])
      vec <- c(vec, num)
    }
    selmatrix <- transinput[vec,]

    # Get mean and sd of new selected gene matrix
    row.mean <- apply(selmatrix, 1, mean)
    row.sd <- apply(selmatrix, 1, stats::sd)

    # Repeat steps and fit second polynomial to data
    fit <- stats::lm(row.sd ~ poly(row.mean, 2, raw=TRUE))
    res <- fit$residuals
    numgenesel <- length(which(res > 0))
    selgenes <- names(which(res > 0))

    # Obtain selected gene matrix of genes above the second polynomial fit
    vec <- NULL
    for (i in 1:length(selgenes)) {
      num <- which(rownames(selmatrix) == selgenes[i])
      vec <- c(vec, num)
    }
    selmatrix2 <- selmatrix[vec,]

    # Obtain mean and sd for new genes above the second polynomial fit
    row.mean <- apply(selmatrix2, 1, mean)
    row.sd <- apply(selmatrix2, 1, stats::sd)

    # Fit a third polynomial to the data
    fit <- stats::lm(row.sd ~ poly(row.mean, 2, raw=TRUE))
    res <- fit$residuals
    numgenesel <- length(which(res > 0))
    selgenes <- names(which(res>0))

    # Obtain selected gene matrix of genes above the third polynomial fit
    vec <- NULL
    for (i in 1:length(selgenes)) {
      num <- which(rownames(selmatrix) == selgenes[i])
      vec <- c(vec, num)
    }

    # Return the number of genes
    gene_poly <- length(vec)
    print(paste("The poly gene probe number is: ", gene_poly, sep=""))
    return(gene_poly)

  }

  # For adaptive method
  if (Adaptive == TRUE) {

    ## Code for Gaussian Mixture Modeling and Selection Process ##
    ## PART 1: NORMALIZATION OF DATA  ##

    # Function to normalize data to bring values into alignment
    nor.min.max <- function(x) {
      x.min <- min(x)
      x.max <- max(x)
      x <- (x - x.min) / (x.max - x.min)
      return (x)
    }

    # Makes Expression Data Numeric
    colname <- colnames(data.exp)
    data.exp <- t(apply(data.exp, 1, as.numeric))
    colnames(data.exp) <- colname

    # Normalization of expression data
    data.min.max <- t(apply(data.exp, 1, nor.min.max))
    rownames(data.min.max) <- rownames(data.exp)
    print("Your data has been normalized")

    ######## PART 2 Generate Simulated Data #########

    # Make normalized data a numeric matrix
    colname <- colnames(data.min.max)
    data.min.max <- t(apply(data.min.max, 1, as.numeric))
    colnames(data.min.max) <- colname

    #Stores number of columns and rows of normalized data into variables
    dsNoOfCol <- NCOL(data.min.max)
    dsNoOfRow <- NROW(data.min.max)

    #specify the names of columns in the data frame
    colnames(data.min.max) <- c(1:dsNoOfCol)
    data.min.max <- data.min.max[, 1:dsNoOfCol]

    #apply function will apply sd and mean to each value in the object
    tmp.mean <- apply(data.min.max, 1, mean)
    tmp.sd <- apply(data.min.max, 1, stats::sd)

    # rnorm generates a random distribution of numbers
    data.sim <- matrix(data=0, nrow=dsNoOfRow, ncol=dsNoOfCol)
    rownames(data.sim) <- rownames(data.min.max)

    for (i in 1:dsNoOfRow) {
      currGeneExpLevel <- stats::rnorm(dsNoOfCol, tmp.mean[i], tmp.sd[i])
      data.sim[i,] <- (currGeneExpLevel)
    }
    print("The simulated data file has been finished")

    ###### PART 3A Apply Gaussian Mixture Modeling to Real Data ######

    # Define file prefix name and other parameters for GMM
    filePrefix <- paste(input, "real_", sep=".")
    expColStart <- 1
    MAX_NO_OF_DIGITS <- 20

    dsNoOfCol <- NCOL(data.min.max)
    dsNoOfRow <- NROW(data.min.max)

    #Set column names as "number 1" through value of dsNoOfCol
    colnames(data.min.max) <- c(expColStart:dsNoOfCol)

    #Stores column numbers into a data frame
    data.min.max <- data.min.max[, expColStart:dsNoOfCol]

    #paste combines the vectors into one string
    #("number of samples: 234" for example)
    print(paste("number of samples: ", dsNoOfCol, sep=""))
    print(paste("number of genes: ", dsNoOfRow, sep=""))

    # calculating the variances
    print(paste("Begin applying Gaussian Mixture Modeling to real data",
                "may take several minutes"))

    for (i in 1:dsNoOfRow) {
      print(i)
      currGeneExpLevel <- data.min.max[i, expColStart:dsNoOfCol]

      # model the gene using Mclust
      gClust <- mclust::Mclust(t(currGeneExpLevel), G=1:6)

      output <- file(paste(filePrefix, "output_pro.txt", sep=""), "at")
      cat(rownames(data.min.max)[i], format(gClust$parameters$pro,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefix, "output_mean.txt", sep=""), "at")
      cat(rownames(data.min.max)[i], format(gClust$parameters$mean,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefix, "output_modelNameG.txt",
          sep=""), "at")
      cat(rownames(data.min.max)[i], gClust$G, gClust$modelName,
          file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefix, "output_sigmasq.txt",
          sep=""), "at")
      cat(rownames(data.min.max)[i],
          format(gClust$parameters$variance$sigmasq,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefix, "output_classification.txt",
          sep=""), "at")
      cat(rownames(data.min.max)[i], gClust$classification,
          file=output, sep="\t")
      cat("\n", file=output)
      close(output)
    }

    print("Gaussian mixture modeling has been applied to your real data")

    ########## PART 3B: Apply GMM to Simu Data  ############

    # Define prefix name for simu data and other parameters for GMM
    expColStart <- 1
    filePrefixsim <- paste(input, "Simu_", sep=".")
    MAX_NO_OF_DIGITS <- 20

    dsNoOfCol <- NCOL(data.sim)
    dsNoOfRow <- NROW(data.sim)

    #Set column names as "number 1" through value of dsNoOfCol
    colnames(data.sim) <- c(expColStart:dsNoOfCol)

    #Stores column numbers into a data frame
    data.sim <- data.sim[, expColStart:dsNoOfCol]

    print(paste("number of samples: ", dsNoOfCol, sep=""))
    print(paste("number of genes: ", dsNoOfRow, sep=""))

    print(paste("Begin applying Gaussian Mixture Modeling to Simulated",
                "Data, may take several minutes"))

    for (i in 1:dsNoOfRow) {
      print(i)
      currGeneExpLevel<- data.sim[i,expColStart:dsNoOfCol]

      # model the gene using Mclust
      gClust <- mclust::Mclust(t(currGeneExpLevel), G=1:6)

      output <- file(paste(filePrefixsim, "output_pro.txt",
                           sep=""), "at")
      cat(rownames(data.sim)[i], format(gClust$parameters$pro,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefixsim, "output_mean.txt",
          sep=""), "at")
      cat(rownames(data.sim)[i], format(gClust$parameters$mean,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefixsim, "output_modelNameG.txt",
          sep=""), "at")
      cat(rownames(data.sim)[i], gClust$G, gClust$modelName,
          file=output, sep="\t")
      cat("\n", file=output)
      close(output)

      output <- file(paste(filePrefixsim, "output_sigmasq.txt",
          sep=""), "at")
      cat(rownames(data.sim)[i],
          format(gClust$parameters$variance$sigmasq,
          digits=MAX_NO_OF_DIGITS), file=output, sep="\t")
      cat("\n", file = output)
      close(output)

      output <- file(paste(filePrefixsim, "output_classification.txt",
          sep=""), "at")
      cat(rownames(data.sim)[i], gClust$classification,
          file=output, sep="\t")
      cat("\n", file = output)
      close(output)
    }

    print(paste("Gaussian Mixture Modeling has been applied to",
                "your simulated data"))

    ####### PART 4: Determine GMM/TRank estimated number of Genes #######

    # Function that sorts a matrix
    mat.sort <- function(mat,n) {
      mat[rank(as.numeric(mat[,n]),
               ties.method="random"),] <- mat[c(1:nrow(mat)),]
      return(mat)
    }

    # Function to compute the TRank value
    compute.TRank <- function(fn.dataModelNameG,
        fn.dataClassification, fn.dataSigmasq, fn.dataMean,
        fn.dataPro,curr.Sigmasq) {
      # Model Name and G
      dataModelNameG <- utils::read.delim2(fn.dataModelNameG,
          header=FALSE, row.names=1, sep="\t", quote="\"",
          dec=".", fill=TRUE)

      # discrete classification produced by mclust
      dataClassification <- utils::read.delim2(fn.dataClassification,
          header=FALSE, row.names=1, sep="\t", quote="\"",
          dec=".", fill=TRUE)

      # Extract Properties of the input file
      # Number of genes
      noOfGene <- NROW(dataModelNameG)
      # Maximum number of gaussians
      maxG <-max(dataModelNameG[, 1])
      # Number of samples
      noOfSample <- NCOL(dataClassification)

      # sigma square produced by mclust
      dataSigmasq <- utils::read.delim2(fn.dataSigmasq,
          header=FALSE, row.names=1, col.names=c(1:(maxG + 1)),
          sep="\t", quote="\"", dec=".", fill=TRUE)

      # mean produced by mclust
      dataMean <- utils::read.delim2(fn.dataMean,
          header=FALSE, row.names=1, col.names=c(1:(maxG + 1)),
          sep="\t", quote="\"", dec=".", fill=TRUE)

      # pi-zero produced by mclust
      dataPro <- utils::read.delim2(fn.dataPro,
          header=FALSE, row.names=1, col.names=c(1:(maxG + 1)),
          sep="\t", quote="\"", dec=".", fill=TRUE)

      # discrete classification produced by mclust
      dataClassification <- utils::read.delim2(fn.dataClassification,
          header=FALSE, row.names=1, col.names=c(1:(noOfSample + 1)),
          sep="\t", quote="\"", dec=".", fill=TRUE)

      #-------------------------------------------------------
      # This matrix stores the gene selection criteria
      #-------------------------------------------------------
      allCriteria <- matrix(nrow=noOfGene, ncol=8)
      print("Begin computing TRank value")
      #-------------------------------------------------------
      #	Calculating the gene selection criteria
      #-------------------------------------------------------
      for (i in 1:noOfGene) {
        print(paste("Now processing Gene ", i))

        # Get the model information.
        # Extract the number of clusters and model name for each gene
        tmpMMADBID <- as.character(rownames(dataModelNameG[i,]))
        tmpG <- dataModelNameG[i, 1]
        tmpModelName <- as.character(dataModelNameG[i, 2])

        # Temporary MATRIX variables to store the statistical criteria
        tmpPWT <- matrix(data=0, nrow=tmpG, ncol=tmpG)
        #-----------------------------------------------------------
        # Calculate the MIN & MAX proportion (pi-zero)
        #-----------------------------------------------------------
        if(tmpG > 1) {
          tmpMinPro <- min(as.numeric(dataPro[tmpMMADBID, 1:tmpG]))
          tmpMaxPro <- max(as.numeric(dataPro[tmpMMADBID, 1:tmpG]))

          #Added to capture the min and max mean (15 July 2011)
          tmpMinMean <- min(as.numeric(dataMean[tmpMMADBID, 1:tmpG]))
          tmpMaxMean <- max(as.numeric(dataMean[tmpMMADBID, 1:tmpG]))
        }

        else {
          tmpMinPro <- 0
          tmpMaxPro <- 0

          #Added to capture the min and max mean (15 July 2011)
          tmpMinMean <- min(as.numeric(dataMean[tmpMMADBID, 1]))
          tmpMaxMean <- max(as.numeric(dataMean[tmpMMADBID, 1]))
        }

        #-----------------------------------------------
        # Experiment (9) - sum(pair-wised)
        #-----------------------------------------------

        for(j in 1:(tmpG - 1)) {
          if(tmpModelName == "X") {
            # do nothing here
          }
          else {
            for(k in (j + 1):(tmpG - 1)) {
              # added on 15 July 2011, otherwise wrong.
              if(k > j) {
                if(tmpModelName == "E") {
                  tmpPWT[j, k] <- (abs(dataMean[tmpMMADBID,j]
                      - dataMean[tmpMMADBID,k]) /
                      sqrt(2 * dataSigmasq[tmpMMADBID, 1] +
                      curr.Sigmasq)) * ((dataPro[tmpMMADBID, j] *
                      dataPro[tmpMMADBID, k]) /
                      ((dataPro[tmpMMADBID, j] +
                      dataPro[tmpMMADBID, k]) ^ 2))
                }
                else if(tmpModelName == "V") {
                  tmpPWT[j, k] <- (abs
                      (dataMean[tmpMMADBID, j] -
                      dataMean[tmpMMADBID, k]) /
                      sqrt(dataSigmasq[tmpMMADBID, j] +
                      dataSigmasq[tmpMMADBID, k] +
                      curr.Sigmasq)) * ((dataPro[tmpMMADBID, j] *
                      dataPro[tmpMMADBID, k]) /
                      ((dataPro[tmpMMADBID, j] +
                      dataPro[tmpMMADBID, k]) ^ 2))
                }
              }
            }
          }
        }

        # store all criteria
        allCriteria[i, 1] <- tmpMMADBID
        allCriteria[i, 2] <- tmpG
        allCriteria[i, 3] <- tmpModelName
        allCriteria[i, 4] <- tmpMinPro
        allCriteria[i, 5] <- tmpMaxPro
        #allCriteria[i, 6] = max(as.numeric(tmpPWT))# max
        allCriteria[i, 6] <- sum(as.numeric(tmpPWT))	# sum
        #Added to capture the min and max mean (15 July 2011)
        allCriteria[i, 7] <- tmpMinMean
        allCriteria[i, 8] <- tmpMaxMean
      }

      #-----------------------------------------------------------
      # sort the matrix by one of the allCriteria
      #-----------------------------------------------------------
      sortedAllCriteria <- mat.sort(allCriteria, 6)
      return (sortedAllCriteria)
    }

    #-------------------------------------------
    #	Global constants
    #-------------------------------------------
    MAX_NO_OF_DIGITS <- 20

    # Read in simu data files for TRank calculation
    tmp.fn.dataModelNameG <- paste(filePrefixsim,
        "output_modelNameG.txt", sep="")
    tmp.fn.dataClassification <- paste(filePrefixsim,
        "output_classification.txt", sep="")
    tmp.fn.dataSigmasq <- paste(filePrefixsim,
        "output_sigmasq.txt", sep="")
    tmp.fn.dataMean <- paste(filePrefixsim,
         "output_mean.txt", sep="")
    tmp.fn.dataPro <- paste(filePrefixsim,
         "output_pro.txt", sep="")

    # Calculate the fudge factor for the simu data
    aS0 <- NULL
    currData <- NULL
    for(i in 1:NROW(data.sim)) {
      currData <- as.numeric(data.sim[i,])
      aS0 <- c(aS0,stats::var(currData))
    }

    tmp.Sigmasq <- stats::median(aS0)
    print(paste("The simulated sigamasq is: ", tmp.Sigmasq, sep=""))

    # Call TRank function to calculate simulated data TRank scores
    x <- compute.TRank(tmp.fn.dataModelNameG, tmp.fn.dataClassification,
        tmp.fn.dataSigmasq, tmp.fn.dataMean, tmp.fn.dataPro, tmp.Sigmasq)
    rownames(x) <- as.character(t(x[, 1]))

    # Begin reading in real data files for TRank calculation
    tmp.fn.dataModelNameG <- paste(filePrefix,
        "output_modelNameG.txt", sep="")
    tmp.fn.dataClassification <- paste(filePrefix,
        "output_classification.txt", sep="")
    tmp.fn.dataSigmasq <- paste(filePrefix,
        "output_sigmasq.txt", sep="")
    tmp.fn.dataMean <- paste(filePrefix,
        "output_mean.txt", sep="")
    tmp.fn.dataPro <- paste(filePrefix,"output_pro.txt", sep = "")

    #calculating the fuge factor s0 for real data
    aS0 <- NULL
    currData <- NULL

    for(i in 1:NROW(data.min.max)) {
      currData <- as.numeric(data.min.max[i,])
      aS0 <- c(aS0,stats::var(currData))
    }

    tmp.Sigmasq <- stats::median(aS0)
    print(paste("The real sigamasq is: ", tmp.Sigmasq, sep=""))

    # Call TRank function to calculate TRank scores for real data
    y <- compute.TRank(tmp.fn.dataModelNameG, tmp.fn.dataClassification,
        tmp.fn.dataSigmasq, tmp.fn.dataMean, tmp.fn.dataPro, tmp.Sigmasq)
    rownames(y) <- as.character(t(y[, 1]))

    # Store TRank scores for simu (x.score) and real data (y.score)
    x.score <- as.numeric(t(x[, 6]))
    y.score <- as.numeric(t(y[, 6]))

    #FDR calculations
    maxScore <- max(c(x.score, y.score))
    currCutoff <- min(c(x.score, y.score))
    FDR <- NULL
    bin.size <- 0.001

    while(currCutoff <= maxScore) {
      a <- sum(x.score > currCutoff)
      b <- sum(y.score > currCutoff)
      currFDR <- a / b
      FDR <- rbind(FDR, c(a, b, currCutoff, currFDR))
      currCutoff <- currCutoff + bin.size
    }

    # Choose the TRank with a FDR value of 0.01 for the lowest number
    # of selected genes a represents the index of the TRank
    # value to be chosen

    FDRvalue <- 0.01
    index <- which(FDR[, 4] >= FDRvalue)

    # Conditional if FDR value never reaches specified value
    # Reduce the FDRvalue

    if (length(index) == 0) {
      index <- which(FDR[, 4] >= (FDRvalue / 10))
      if (length(index) == 0) {
        index <- which(FDR[, 4] >= (FDRvalue / 50))
        if (length(index) == 0) {
          index <- which(FDR[, 4] >= (FDRvalue / 100))
          if (length(index) == 0) {
            index <- which(FDR[, 4] >= (FDRvalue / 1000))
          }
        }
      }
    }

    genenum <- max(index)
    TRank <- FDR[genenum, 3]
    # Conditionals to make sure only one TRank value is chosen
    if (length(TRank) == 0) {
      TRank <- 0.30
      print("Arbitrary TRank value of 0.30 chosen")
      print("The specified FDR cutoff did not occur in the data")
    }

    # Conditional if more than 1 TRank value exists with the FDR cutoff
    if (length(TRank) > 1) {
      TRank <- min(TRank)
    }
    # Save the TRank determined genes to object
    sel.Probe <- as.character(y[as.numeric(y[, 6]) >= TRank, 1])
    gene_num3 <- length(sel.Probe)
    print(paste("The adaptive gene probe number is: ", gene_num3, sep=""))
    return(gene_num3)
  }

}

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)

test_number_probes <- function() {
    checkException(number_probes(input="test", data.exp="string", Fixed=1,
        Percent=NULL, Poly=NULL, Adaptive=NULL),
        msg="Unable to use string as data.exp")
    checkException(number_probes(input="test", data.exp=test_matrix, Fixed=NULL,
        Percent=TRUE, Poly=TRUE, Adaptive=NULL),
        msg="Can't have multiple options as true")
}
