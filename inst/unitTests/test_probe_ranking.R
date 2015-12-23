probe_ranking <- function(input, probe_number,
    probe_num_selection="Fixed_Probe_Num", data.exp, method="SD_Rank") {

    WriteMatrixToFile <- function(tmpMatrix, tmpFileName,
        blnRowNames, blnColNames) {
    output <- file(tmpFileName, "at")
    utils::write.table(tmpMatrix, output, sep="\t", quote=FALSE,
        row.names=blnRowNames, col.names=blnColNames)
    close(output)
    }

  #Makes Expression Data Numeric
  colname <- colnames(data.exp)
  data.exp <- t(apply(data.exp, 1, as.numeric))
  colnames(data.exp) <- colname

  # If CV_Rank option is chosen
  if (method == "CV_Rank") {
    #Mean and SD Filter
    row.sd <- apply(data.exp, 1, stats::sd)
    row.mean <- apply(data.exp, 1, mean)
    row.cv <- row.sd / row.mean

    # Obtain new gene expression matrix of CV-ranked genes
    data.exp <- data.frame(data.exp, row.cv)
    data.exp <- data.exp[order(data.exp[, dim(data.exp)[2]],
        decreasing=TRUE),]
    data.exp <- data.exp[, -dim(data.exp)[2]]
    newexp.cv <- data.exp[1:probe_number,]
    WriteMatrixToFile(newexp.cv, paste(input, probe_num_selection,
        "CV_Rank_Selected_Gene_Exp.txt", sep="."), TRUE, TRUE)
    print(paste("The selected CV_Rank Gene Expression text file",
        "has been written"))
    return(newexp.cv)
  }

  # If CV_Guided method is chosen
  if (method == "CV_Guided") {
    # Shift the data by the minimum value
    num <- min(data.exp)
    transinput <- data.exp - num

    # Obtain selected gene matrix
    constant <- stats::sd(transinput) / mean(transinput)
    row.sd <- apply(transinput, 1, stats::sd)
    row.mean <- apply(transinput, 1, mean)
    x <- 0
    y <- 55000
    while (y > probe_number) {
      y <- length(which(row.mean > (0.0001 *x) &
          row.sd > (0.0001* constant *x)))
      x <- x + 1
    }
    x <- x - 1

    selmatrix <- transinput[row.mean > (0.0001 *x) &
          row.sd > (.0001 * constant * x),]
    WriteMatrixToFile(selmatrix, paste(input, probe_num_selection,
          "CV_Guided_Selected_Gene_Exp.txt", sep='.'), TRUE, TRUE)
    print(paste("The selected CV_Guided Gene Expression",
          "text file has been written"))
    return(selmatrix)
  }

  if (method == "SD_Rank") {
    # Shift the data by the minimum value
    num <- min(data.exp)
    transinput <- data.exp - num

    # Get selected gene expression matrix
    row.sd <- apply(transinput, 1, stats::sd)
    transinput <- data.frame(transinput, row.sd)
    transinput <- transinput[order(transinput[, dim(transinput)[2]],
        decreasing = TRUE),]
    transinput <- transinput[, -dim(transinput)[2]]
    newexp.sd <- transinput[1:probe_number,]

    WriteMatrixToFile(newexp.sd, paste(input, probe_num_selection,
        "SD_Rank_Selected_Gene_Exp.txt", sep="."), TRUE, TRUE)
    print(paste("The selected SD_Rank Gene Expression",
        "text file has been written"))
    return(newexp.sd)
  }

  if (method == "Poly") {
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
    for (i in 1:length(selgenes)) {
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

    # Order the genes by highest res number
    df <- data.frame(names(res), res)
    newdf <- df[order(df[,dim(df)[2]], decreasing=TRUE),]

    # Conditional if the user input gene number is larger than the
    # Poly selected amount produce an error prompt
    if (dim(newdf)[1] < probe_number) {
      print(paste("The number of genes determined by the Poly",
          "filter was: ", dim(newdf)[1], sep=""))
      stop("Please choose a gene number equal to or less
          than the determined gene number")
    }

    poly_genes <- newdf[1:probe_number,]
    # Obtain the selected gene matrix for genes above the third polynomial
    selgenes <- rownames(poly_genes)

    vec <- NULL
    for (i in 1:length(selgenes)) {
      num <- which(rownames(selmatrix2) == selgenes[i])
      vec <- c(vec, num)
    }
    selmatrix3 <- selmatrix2[vec,]
    WriteMatrixToFile(selmatrix3, paste(input,
        "Poly_Selected_Gene_Exp.txt", sep="."), TRUE, TRUE)
    print(paste("The selected Poly Gene Expression text",
        "file has been written"))
    return(selmatrix3)
    }

  # Conditional to produce error if user does not choose one
  # of available options
  choices <- c("CV_Rank", "CV_Guided", "SD_Rank", "Poly")
  if (is.na(match(method, choices)) == TRUE) {
    stop("Please choose one of the available methods: CV_Rank,
         CV_Guided, SD_Rank, or Poly")
  }
}

# Make a test matrix
test_data <- c(1:8)
test_matrix <- matrix(data=test_data, nrow=2, ncol=4)

test_probe_ranking <- function() {
    checkException(probe_ranking(input="string", probe_number="word",
        probe_num_selection="Fixed_Probe_Num", data.exp=test_matrix,
        method="SD_Rank"), msg="Unable to use string as probe_number")
}
