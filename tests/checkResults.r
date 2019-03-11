
print2 = function(string){
    cat(string)
    cat("\n")
}

trueGenotypes = as.matrix(read.table("baseData/trueGenotypes.txt"))

getGwasCor = function(mat, true) {
    cors = sapply(1:ncol(true), function(ii){
        cor(true[,ii], mat[,ii], use = "pair")
    })
    return(mean(cors, na.rm = T))
}

assessPeeling = function(filePrefix){

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))
    oldFile = NULL
    tryCatch({
        oldFile = as.matrix(read.table(paste0("lastStable/", filePrefix)))
        }, error = function(err){}, warning=function(err){})
    if(is.null(oldFile)) {
        print("NO OLD FILE")
        oldFile = newFile
    }


    print2(" ")
    print2(paste("Assessing peeling file:", filePrefix))
    print2(paste("Checking if outputs are equal:", all(newFile == oldFile)))

    newAcc = getGwasCor(newFile[,-1], trueGenotypes[,-1])
    oldAcc = getGwasCor(oldFile[,-1], trueGenotypes[,-1])
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))

}

assessPeeling("peeling.dosages")
assessPeeling("peeling.multi.seq.dosages")
assessPeeling("peeling.single.dosages")
assessPeeling("peeling.hybrid.dosages")
assessPeeling("peeling.hybrid.seq.dosages")


pedigree = read.table("baseData/pedigree.txt")
assessAssign = function(filePrefix){
    lastGen = 801:1000

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))
    oldFile = NULL
    tryCatch({
        oldFile = as.matrix(read.table(paste0("lastStable/", filePrefix)))
        }, error = function(err){}, warning=function(err){})
    if(is.null(oldFile)) {
        print("NO OLD FILE")
        oldFile = newFile
    }


    print2(" ")
    print2(paste("Assessing Assign file:", filePrefix))
    print2(paste("Checking if outputs are equal:", all(newFile == oldFile)))

    newAcc = mean(newFile[lastGen,2] == pedigree[lastGen,2])
    oldAcc = mean(oldFile[lastGen,2] == pedigree[lastGen,2])
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))

}
assessAssign("assign.pedigree")
assessAssign("assign.noDam.pedigree")
assessAssign("assign.seq.pedigree")


pedigreeExtended = read.table("baseData/pedigree.extended")
assessMGS = function(filePrefix){
    lastGen = 801:1000

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))
    oldFile = NULL
    tryCatch({
        oldFile = as.matrix(read.table(paste0("lastStable/", filePrefix)))
        }, error = function(err){}, warning=function(err){})
    if(is.null(oldFile)) {
        print("NO OLD FILE")
        oldFile = newFile
    }

    print2(" ")
    print2(paste("Assessing MGS file:", filePrefix))
    print2(paste("Checking if outputs are equal:", all(newFile == oldFile)))

    newAcc = mean(newFile[,2] == pedigreeExtended[lastGen,4])
    oldAcc = mean(oldFile[,2] == pedigreeExtended[lastGen,4])
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))

}
assessMGS("mgsAssign.grandsires")



assessImputation = function(filePrefix){

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))

    oldFile = NULL
    tryCatch({
        oldFile = as.matrix(read.table(paste0("lastStable/", filePrefix)))
        }, error = function(err){}, warning=function(err){})
    if(is.null(oldFile)) {
        print("NO OLD FILE")
        oldFile = newFile
    }

    print2(paste("Assessing imputation file:", filePrefix))
    print2(paste("Checking if outputs are equal:", all(newFile == oldFile, na.rm=T)))


    newFile[newFile == 9] = NA
    oldFile[oldFile == 9] = NA

    oldCor = getGwasCor(oldFile[,-1], trueGenotypes[,-1])
    oldAcc = mean(oldFile[,-1] == trueGenotypes[,-1], na.rm=T)
    oldYield = mean(!is.na(oldFile[,-1]))

    newCor = getGwasCor(newFile[,-1], trueGenotypes[,-1])
    newAcc = mean(newFile[,-1] == trueGenotypes[,-1], na.rm=T)
    newYield = mean(!is.na(newFile[,-1]))

    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))
    print2(paste("Comparing yield: ", round(oldYield, digits=3), "->", round(newYield, digits=3)))
    print2(paste("Comparing correlations: ", round(oldCor, digits=3), "->", round(newCor, digits=3)))

}
# assessImputation("genotypes.txt")
assessImputation("imputation.genotypes")






