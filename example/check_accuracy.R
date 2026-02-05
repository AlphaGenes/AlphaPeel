#!/usr/bin/env Rscript

trueGenotypes = as.matrix(read.table("data/trueGenotypes.txt"))

getMarkerCorrelations = function(mat, true) {
  cors = sapply(1:ncol(true), function(ii) {
    cor(true[, ii], mat[, ii], use = "pair")
  })
  return(mean(cors, na.rm = TRUE))
}

assessPeeling = function(filePrefix) {
  newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))
  print(" ")
  print(paste("Assessing peeling file:", filePrefix))
  newAcc = getMarkerCorrelations(newFile[, -1], trueGenotypes[, -1])
  print(paste("Comparing accuracies: ", round(newAcc, digits = 3)))
}

assessPeeling("multilocus.dosage.txt")
assessPeeling("hybrid.dosage.txt")
