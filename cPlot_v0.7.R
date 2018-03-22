# cPlot2.R
# Thomas Luh
# Sets up functions for a dimensional analysis of cell differentiation, 
# pre-processing of data, t-SNE analysis, and then assignment of all cells
# to a cell lineage. Based on R scripts by Wenji Ma.
# 
# Workflow:
# 0) set working directory to proper location, and load this file into R
# 1) add data via addCB(), addPatient(), addPreNorm()
# 2) (optional) writeOut() the processed/normalized data
# 3) generateOneMap() to check general shape and adjust perplexity/theta
# 4) generateManyMaps() and hope for a map with good orientation
# 4.5) manually apply rotation/flip transformations to t-SNE coordinates
# 5) pickMap() to import a file with desirable t-SNE coordinates
# 6) assignLineages() to all sample cells
# 7) compileAllData() to make an easily-readable compilation of all data

# version history
# v0.7: public release

readData <- function (inFile, keyword = "") {
# This function accepts a .csv as input, and the creates a dataframe with 
# named rows.
# 
# args:
# inFile: Filepath to .csv file.
# keyword: String of characters prepended to each row's name.
# 
# returns:
# Dataframe with named rows.

  data <- read.csv(inFile, header = TRUE, check.names = TRUE)
  rownames(data) <- paste0(keyword, "_", rownames(data), "-", data$Progenitor)
  return(data)
}

normalizeData <- function(counts, method = NULL) {
# This function normalizes a dataframe according to a protocol of choice.
# cPlot2 only uses the DESeq method.
# DESeq method requires that each row is named, and that the progenitor name
# is present in the row names!
# 
# args:
# counts: The data to be normalized.
# method: The normalization method to be used. Options are "median", "edgeR", 
#         "DESeq", and "sum".
# 
# returns:
# Normalized dataframe.

  if (method=="median"){
    normalizedCounts <- apply(counts, 2, function(x){x / median(x[x > 5])})
  }
  if (method=="edgeR"){
    library(edgeR)
    normalizedCounts <- t(t(counts) * calcNormFactors(counts))
  }
  if (method=="DESeq"){
    library(DESeq2)
    normalizedCounts <- t(t(counts) / estimateSizeFactorsForMatrix(counts))
  }
  if (method=="sum"){
    normalizedCounts <- apply(counts, 2, function(x){x / sum(x)})
  }
  return(normalizedCounts)
}

normReadoutByPatient <- function(dataIn, dataCountsIn) {

  dataCountsIn_t <- t(dataCountsIn[, c(-1, -2)])
  dataCountsIn_t <- log10(as.matrix(dataCountsIn_t) + 1)
  dataCountsIn_norm <- dataCountsIn_t
  
  progenitorTypes <- unique(dataCountsIn$Progenitor)
  
  # for each type of progenitor:
  for(progenitor in progenitorTypes){
    # find mother cells of that progenitor type,
    cell.index <- grepl(paste0("-", progenitor), rownames(dataCountsIn))
    # and then get each mother's batch ID and progeny cell counts.
    batch.vec <- dataIn$Batch[cell.index]
    data.input <- dataCountsIn_t[, cell.index]
    # If the progenitor comes from 2+ batches, then normalize:
    if(length(unique(batch.vec)) <= 1){
      next
    } else {
    dataCountsIn_norm[, cell.index] <- 
      normBatchGeometricMean(dataCountsIn_t[, cell.index], batch.vec)
    }
  }
  
  dataCountsOut_norm <- 
    data.frame(Progenitor = dataIn$Progenitor, t(dataCountsIn_norm))
  dataCountsOut <- 
    data.frame(Progenitor = dataIn$Progenitor, 10^(dataCountsOut_norm[, -1]) -1)
  return(dataCountsOut)  
}

normBatchGeometricMean <- function(data.input, batch.vec){
  #
  #
  # args:
  # data.input: a log10-transformed matrix; each column = one sample
  # batch.vec:  a vector containing batch IDs of samples
  #
  # returns:
  # a modified dataset that has been normalized according to geometric means
  
  # remove duplicate subjects from batch.vec
  batches <- unique(batch.vec)
  # set up an empty statistics.matrix
  statistics.matrix <- 
    matrix(0, nrow = nrow(data.input), ncol = length(batches))
  
  # for each type of progeny and each batch:
  for(i in 1:nrow(data.input)){
    for(j in 1:length(batches)){
      # calculate the mean of the cell counts
      statistics.matrix[i, j] <- mean(data.input[i, batch.vec == batches[j]])
    }
  }
  
  # calculate the difference from each reading to the median per progeny
  ref.sample <- apply(statistics.matrix, 1, median)
  ratio.matrix <- statistics.matrix - ref.sample
  # then, find the batches' size factors (median difference for each batch)
  size.factor <- apply(ratio.matrix, 2, median)
  
  # save the size factors
  size.factor.vec <- rep(0, length(batch.vec))
  for(i in 1:length(batches)){
    size.factor.vec[batch.vec == batches[i]] <- size.factor[i]
  }
  # subtract the size factor from each sample
  
  out <- data.input - size.factor.vec
  return(out)
}

processData <- function(dataIn, firstData, dataCols,
                        normCountsBySubject = TRUE,
                        transformData = TRUE,
                        normProgenies = TRUE) {
  # Here, we process and normalize data. This is the "master" normalization
  # function.
  #
  # args:
  # dataIn: the name of a data.frame to be processed.
  # firstData: the first column containing data. First column = 1.
  # dataCols: the number of columns containing data.
  # normProgenies: ??
  # normCountsBySubject: ??
  # transformData: remove low cell counts and compensate for cell loss.
  # 
  # returns:
  # processed dataframe, containing ONLY columns of data (i.e. non-data 
  # columns are purged from the dataframe)
  
  lastData <- firstData + dataCols - 1      # the last column containing data
  
  # create the local dataframe.
  dataOut <- dataIn[, c(match(c("Progenitor","Batch"), names(dataIn)), 
                        (firstData):(lastData))]
  rownames(dataOut) <- rownames(dataIn)
  colnames(dataOut) <- c("Progenitor", 
                         "Batch", 
                         colnames(dataIn)[(firstData):(lastData)])
  
  if(normCountsBySubject){
    dataOut <- normReadoutByPatient(dataIn, dataOut)
  }
  
  if(transformData){
    
    # for CDP cells, neglect low cell counts
    CDPcutoff <- 2
    for(i in 2:ncol(dataOut)){
      dataOut[which(dataOut[, i] <= CDPcutoff), i] <- 0
    }
    
    # for non-CDP cells, neglect low cell counts
    nonCDPcutoff <- 6
    notCDP <- !grepl("CDP", dataOut[, 1])
    for(i in 2:ncol(dataOut)){
      dataOut[((notCDP) & (dataOut[, i] <= nonCDPcutoff)), i] <- 0
    }
    
    # increase cell count for all cells
    cellLossFactor <- 4.5
    dataOut[, -1] <- dataOut[, -1] * cellLossFactor
  }
  
  if(normProgenies){
    balanceProgeny <- normalizeData(dataOut[, -1], method = "DESeq")
    dataOut <- data.frame(Progenitor = dataOut$Progenitor, balanceProgeny)
  }
  
  return(dataOut)
}

addCB <- function (inFile, firstData, dataCols, keyword = "") {
  # This function processes data from cord blood sources. Make sure that
  # the input file has the same format (especially column names!)
  # as already added data.
  #
  # args:
  # inFile: the name of a csv to be pre-processed.
  # firstData: the first column of data. (1-indexed)
  # dataCols: the number of columns of data. (1-indexed)
  # keyword: String of characters prepended to each row's name.
  #
  # returns:
  # Nothing. This is a mutator method.
  
  dataCB <- readData(inFile, keyword)             # read raw data
  dataCountsCB <- processData(dataCB,             # process data
                              firstData = firstData,
                              dataCols = dataCols)
  if(exists("data") && is.data.frame(get("data"))) {    # save raw data
    data <<- rbind(data, dataCB)
  } else {
    assign("data", dataCB, envir=.GlobalEnv)
  }
  if(exists("dataCounts") && is.data.frame(get("dataCounts"))) {
    dataCounts <<- rbind(dataCounts, dataCountsCB)
  } else {
    assign("dataCounts", dataCountsCB, envir=.GlobalEnv)
  }
}

addPatient <- function(inFile, firstData, dataCols, keyword = "") {
  # This function processes data from non-cord-blood sources. Make sure that
  # the input file has the same format (especially column names!)
  # as already added data.
  # args:
  # inFile: the name of a csv to be pre-processed.
  # firstData: the first column of data. (1-indexed)
  # dataCols: the number of columns of data. (1-indexed)
  # keyword: String of characters prepended to each row's name.
  #
  # returns:
  # Nothing. This is a mutator method.
  
  dataPt <- readData(inFile, keyword)             # read raw data
  dataCountsPt <- processData(dataPt,             # process data
                              normProgenies = FALSE,
                              firstData = firstData,
                              dataCols = dataCols)
  # save raw data                                                  
  if(exists("data") && is.data.frame(get("data"))) {    
    data <<- rbind(data, dataPt)
  } else {
    assign("data", dataPt, envir=.GlobalEnv)
  }
  # save processed data
  if(exists("dataCounts") && is.data.frame(get("dataCounts"))) {
    dataCounts <<- rbind(dataCounts, dataCountsPt)  
  } else {
    assign("dataCounts", dataCountsPt, envir=.GlobalEnv)
  }
}

addPreNorm <- function(inFile, firstData, dataCols, keyword = "") {
  # This function processes data from any source, as long as it's been 
  # normalized already. Make sure that the input file has the same format 
  # (especially column names!) as already added data.
  #
  # args:
  # inFile: the name of a csv to be pre-processed.
  # firstData: the first column of data. (1-indexed)
  # dataCols: the number of columns of data. (1-indexed)
  # keyword: String of characters prepended to each row's name.
  #
  # returns:
  # Nothing. This is a mutator method.
  
  lastData <- firstData + dataCols - 1      # the last column containing data
  
  dataNorm <- readData(inFile, keyword)
  dataNorm <- dataNorm[, c(match("Progenitor", names(dataNorm)),
                           (firstData):(lastData))]
  # save raw data                                                  
  if(exists("data") && is.data.frame(get("data"))) {    
    data <<- rbind(data, dataNorm)
  } else {
    assign("data", dataNorm, envir=.GlobalEnv)
  }
  # save processed data
  if(exists("dataCounts") && is.data.frame(get("dataCounts"))) {
    dataCounts <<- rbind(dataCounts, dataNorm)
  } else {
    assign("dataCounts", dataNorm, envir=.GlobalEnv)
  }
}

purify <- function() {
  # This function removes progenitors that output 0 cells from dataCounts.
  # May not work with all data sets.
  # 
  # args:
  # inData: a dataframe to be manipulated
  # 
  # returns:
  # Nothing. This mutates the dataframe.
  dataCounts <- 
    dataCounts[as.logical(dataCounts$G + dataCounts$M + dataCounts$C + 
                            dataCounts$A + dataCounts$P + dataCounts$L),]
  data <- 
    data[as.logical(data$G + data$M + data$C + 
                      data$A + data$P + data$L),]
}

writeOut <- function(outFile = "") {
  # This function exports the currently processed data from all sources.
  #
  # args:
  # outFile: a string to prepend to the filename
  # 
  # returns:
  # Nothing. This creates a .csv in the current working directory.
  
  write.csv(dataCounts, paste0(outFile, "processedData", format(Sys.time(), 
                               "%Y%b%d_%I%M%S%p"), ".csv"))
}

generateOneMap <- function (tsneIn = dataCounts, perIn = 30, thetaIn = 0.1, 
                            randSeed = 42, firstData = 2, dataCols = 6,
                            sampleSize = nrow(dataCounts),
                            iterIn = 1000) {
  # This function generates one plot, based on the user's inputs.
  # 
  # args:
  # tsneIn: dataframe to be manipulated.
  # perIn: perplexity, passed to the Rtsne function.
  # thetaIn: theta, passed to the Rtsne function.
  # randSeed: prepares a random number seed for reproducibility.
  #
  # returns:
  # nothing. generates one plot in current working directory.
  
  library(Rtsne)
  library(ggplot2)
  set.seed(randSeed)
  
  tsneSub <- tsneIn[sample(nrow(tsneIn), sampleSize), ] #random sample
  
  # perform t-SNE analysis
  tsneOut <- Rtsne(as.matrix(tsneSub[,(firstData):(firstData + dataCols - 1)]),
                   perplexity = perIn,
                   theta = thetaIn,
                   max_iter = iterIn,
                   verbose = TRUE,
                   check_duplicates = FALSE)

  # generate and save plot
  mainOut <- paste0("p = ", perIn, ", t = ", thetaIn, ", s = ", randSeed)
  png(paste0("plotSolo","_p", perIn, "_t", thetaIn, "_s", randSeed, 
             "_i", iterIn, ".png"))
  print(qplot(tsneOut$Y[, 1], tsneOut$Y[, 2], main = mainOut))
  # [colour = dataCounts$Progenitor] is a possible color option above
  dev.off()
  
  # generate and save data points
  write.csv(tsneOut$Y, file = paste0("dataSolo", 
                                     "_p", perIn, 
                                     "_t", thetaIn, 
                                     "_s", randSeed, 
                                     "_i", iterIn, ".csv"))
}

generateManyMaps <- function(tsneIn = dataCounts, 
                             firstData = 2, dataCols = 6,
                             perVector = c(50), 
                             thetaVector = c(0.10),
                             iterVector = c(1000),
                             etaVector = c(200),
                             n = 1,
                             seedSet = FALSE,
                             randSeed = 6, 
                             testMode = FALSE,
                             exportAs = "png") {
  # This function generates many plots, and exports them to the current working
  # directory. It also exports the coordinates of the points.
  # 
  # args:
  # tsneIn: dataframe to be manipulated.
  # perIn: perplexity, passed to the Rtsne function.
  # thetaIn: theta, passed to the Rtsne function.
  # n: number of replications of each dataset to be done
  # randSeed: prepares a random number seed for reproducibility.
  # testMode: if on, TSNE runs for 10 iterations (faster).
  # exportAs: select an export method. Options: "png", "pdf" (experimental)
  #
  # returns:
  # Nothing. This creates a series of plots and file outputs in the current 
  # working directory.
  
  library(Rtsne)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  
  set.seed(randSeed)                        # allows reproducibility
  numPlots <- length(iterVector) * 
                length(perVector) * 
                length(thetaVector) * 
                length(etaVector) *
                n                           # number of plots
  curPlot <- 0                              # current plot counter
  lastData <- firstData + dataCols - 1
  
  # reduce calculation time if test mode is active
  if(testMode){iterIn <- 10} else {iterIn <- 1000}
  
  # make the plots
  curPlot <- 0
  cat("Progress: 0 of", numPlots, "\n")
  for(perIn in perVector){
    for(thetaIn in thetaVector){
      for(iterIn in iterVector){
        for(etaIn in etaVector){
            for (i in 1:n){
                if(seedSet){set.seed(randSeed)}      #seed reset for testing
                curPlot <- curPlot + 1                      #increment counter
                tsneOut <- Rtsne(as.matrix(tsneIn[, firstData:lastData]), 
                                 theta = thetaIn,
                                 perplexity = perIn,
                                 max_iter = iterIn,
                                 eta = etaIn,
                                 check_duplicates = FALSE)  # make map
                # generate plots
                mainOut <- paste0("p = ", perIn, ", t = ", thetaIn)
                if(exportAs == "png"){
                    png(paste0("plot", curPlot, "_p", perIn, "_t", thetaIn,
                               "_s", randSeed, "_i", iterIn, ".png"))
                    
                    print(qplot(tsneOut$Y[, 1], tsneOut$Y[, 2], main = mainOut))
                    # colour = dataCounts$Progenitor
                    dev.off()
                } else if (exportAs == "pdf"){
                    #newPlot <- qplot(tsneOut$Y[, 1], tsneOut$Y[, 2], 
                    #                 main = mainOut, 
                    #                 colour = dataCounts$Progenitor)
                    #assign(paste0("plot", curPlot), newPlot)
                }
                
                # generate csv
                write.csv(tsneOut$Y, file = paste0("data", curPlot,"_p", perIn, 
                                                   "_t", thetaIn, "_s", randSeed, 
                                                   "_i", iterIn, ".csv"))
                
                # print status update
                cat("\rProgress:", curPlot, "of", numPlots, "plotted\n")
            }
        }
      }
    }
  }
  
  if(exportAs == "pdf"){
    # second, extend the plot series by increasing counter to a multiple of 15
    # while(nameCount %% 15 != 0){
    # nameCount <- nameCount + 1
    # newPlot <- grid.rect(gp = gpar(col="white"))
    # assign(paste0("plot", nameCount), newPlot)
    #}

    # finally, arrange the plots and generate a .pdf:
    # (this part of the code is broken)
    # pl <- lapply(1:nameCount, function(.x) paste0("plot", .x))
    # ml <- marrangeGrob(pl, nrow=3, ncol=5)
    # ggsave("manyPlotsOutput.pdf", ml, width = 11, height = 8.5)
  }
  
}

pickMap <- function (inFile) {
  # This function lets the user pick a desired t-SNE map, and import it to
  # the workspace.
  # 
  # args:
  # inFile: Filepath to .csv file, written by this script.
  # 
  # returns:
  # A formatted and named dataframe with t-SNE results.
  
  tsneMap <- read.csv(inFile, header = T)
  tsneMap <- tsneMap[2:3]
  rownames(tsneMap) <- rownames(data)
  assign("tsneMap", tsneMap, envir = .GlobalEnv)
}

flipX <- function () {
  tsneMap <- cbind(tsneMap[1], tsneMap[2]*-1)
}

flipY <- function () {
  tsneMap <- cbind(tsneMap[1]*-1, tsneMap[2])
}

rotateClockwise90 <- function () {
  tsneMap <- cbind(tsneMap[2], tsneMap[1]*-1)
}

saveMapb <- function () {
  write.csv(tsneMap, file = paste0(outFile, "rotatedData", format(Sys.time(), 
                               "%Y%b%d_%I%M%S%p"), ".csv"))
}

assignLineages <- function(inData = tsneMap, outFile = "") {
  # This function assigns a lineage to each cell.
  # You must first use pickMap to pick a map, before running this analysis.
  # You must have information about your dataset in the variable dataCounts.
  # This should be leftover from "add" series of functions.
  #
  # args:
  # inData: 
  # outFile: 
  #
  # returns:
  # Nothing. Generates one .csv file in the current working directory.
  
  fractionCutoff <- 0.7
  cellLabels <- colnames(dataCounts)[2:7]
  commitmentMatrix <- dataCounts[, -1] / (rowSums(dataCounts[, 2:7]))
  uniMatrix <- commitmentMatrix > fractionCutoff
  multiPotentProgenitor <- rowSums(commitmentMatrix > fractionCutoff)
  multiPotentProgenitor[is.na(multiPotentProgenitor)] <- 0 # 0 output cells fix
  
  labels <- rep(0, nrow(inData))
  names(labels) <- rownames(inData)
  
  ratio <- labels
  distanceMatrix <- matrix(0, nrow = nrow(inData), ncol = 6)
  rownames(distanceMatrix) <- rownames(inData)
  
  labels <- labels[!is.na(labels)]
  # for unipotent cells, identify the outputs they make the most
  labels[multiPotentProgenitor == 1] <- 
    unlist(apply(commitmentMatrix, 1, function(x){ which(x > fractionCutoff)})) 
  
  fileName <- paste0(outFile, "assignedLineages_FC", fractionCutoff, ".csv")
  if(!file.exists(fileName)){
    for(i in names(multiPotentProgenitor)){
      tempMatrix <- rbind(tsneMap[i, ], tsneMap[multiPotentProgenitor == 1, ])
      distance <- as.matrix(dist(tempMatrix))[1, ]
      
      for(j in 1:6){
        distanceMatrix[i, j] <- 
          min(distance[names(labels)[labels == j & multiPotentProgenitor == 1]],
            na.rm = T)
      }
      labels[i] <- which.min(distanceMatrix[i, ])
      out <- data.frame(lineageAssignment = labels, distanceMatrix)
      colnames(out) <- c("lineageAssignment", paste0("dist_", cellLabels))
      write.csv(out, fileName)
    }
  } else {
    result <- read.csv(fileName, row.names = 1)
    labels <- result$labels
  }
}

compileAllData <- function(inFile, outFile = "") {
  # This function compiles data from (a) the workspace and (b) the file 
  # containing lineage outputs.
  # 
  # args:
  # inFile: Filepath to a clustering result file from assignLineages().
  # outFile: String to prepend to output file.
  # 
  # returns:
  # Nothing. Generates one .csv file in current working directory.

  out <- read.csv(inFile, header = T)[, -1]
  linAssgn <- out$lineageAssignment
  commitment <- rep(0, nrow(dataCounts))
  lineageOut <- commitment
  
  for(i in 1: length(commitment)){
    index <- labels[i] + 1
    commitment[i] <- dataCounts[i, index] / sixLineage[i]
    lineageOut[i] <- dataCounts[i, index]
  }
  
  export <- data.frame(data, dataCounts[, 1:7], tsneMap,  out, lineageOut, 
                       commitment)
  write.csv(export, paste0(outFile, "tsneMapLineageAssignment.csv"), 
            row.names=F)
}