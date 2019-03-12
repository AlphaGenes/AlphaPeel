
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
    print2(" ")
    print2(paste("Assessing peeling file:", filePrefix))
    newAcc = getGwasCor(newFile[,-1], trueGenotypes[,-1])
    print2(paste("Accuracy: ", round(newAcc, digits=3)))

}

assessPeeling("peeling.dosages")
assessPeeling("peeling.multi.seq.dosages")
assessPeeling("peeling.single.dosages")
assessPeeling("peeling.hybrid.dosages")
assessPeeling("peeling.hybrid.seq.dosages")


