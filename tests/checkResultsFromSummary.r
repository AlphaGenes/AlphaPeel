
args = commandArgs(trailingOnly = TRUE)

reference = data.frame(program=c(), file=c(), value = c())
try({
    reference = read.table(args[1], header=T)
})
output = args[2]

print2 = function(string){
    cat(string)
    cat("\n")
}

results = list()

addResult = function(programName, fileName, value){
    code = paste(programName, fileName)
    results[[code]] <<- data.frame(program = programName, file = fileName, value = value)
}
getResult = function(programName, fileName, defaultValue){
    value = reference[reference$program == programName & reference$file == fileName, "value"] 
    if(length(value) == 0){
        # print2(paste("Reference file", programName, fileName, "not found. Returning NA"))
        return(NA)
    }
    return(value)
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
    oldAcc = getResult("Peeling", filePrefix)
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))

    addResult("Peeling", filePrefix, round(newAcc, digits=3))

}

assessPeeling("peeling.dosages")
assessPeeling("peeling.estErrorsAndMaf.dosages")
assessPeeling("peeling.multi.seq.dosages")
assessPeeling("peeling.single.dosages")
assessPeeling("peeling.hybrid.dosages")
assessPeeling("peeling.hybrid.seq.dosages")


pedigree = read.table("baseData/pedigree.txt")
assessAssign = function(filePrefix){
    lastGen = 801:1000

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))

    print2(" ")
    print2(paste("Assessing Assign file:", filePrefix))

    newAcc = mean(newFile[lastGen,2] == pedigree[lastGen,2])
    oldAcc = getResult("Assign", filePrefix)
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))

    addResult("Assign", filePrefix, round(newAcc, digits=3))

}
assessAssign("assign.pedigree")
assessAssign("assign.noDam.pedigree")
assessAssign("assign.seq.pedigree")


pedigreeExtended = read.table("baseData/pedigree.extended")
assessMGS = function(filePrefix){
    lastGen = 801:1000

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))

    print2(" ")
    print2(paste("Assessing MGS file:", filePrefix))

    newAcc = mean(newFile[,2] == pedigreeExtended[lastGen,4])
    oldAcc = getResult("MGS", filePrefix)
    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))
    addResult("MGS", filePrefix, round(newAcc, digits=3))

}
assessMGS("mgsAssign.grandsires")



assessImputation = function(filePrefix){

    newFile = as.matrix(read.table(paste0("outputs/", filePrefix)))

    print2(paste("Assessing imputation file:", filePrefix))


    newFile[newFile == 9] = NA

    oldCor = getResult("ImputeCor", filePrefix)
    oldAcc = getResult("ImputeAcc", filePrefix)
    oldYield =  getResult("ImputeYield", filePrefix)

    newCor = getGwasCor(newFile[,-1], trueGenotypes[,-1])
    newAcc = mean(newFile[,-1] == trueGenotypes[,-1], na.rm=T)
    newYield = mean(!is.na(newFile[,-1]))

    print2(paste("Comparing accuracies: ", round(oldAcc, digits=3), "->", round(newAcc, digits=3)))
    print2(paste("Comparing yield: ", round(oldYield, digits=3), "->", round(newYield, digits=3)))
    print2(paste("Comparing correlations: ", round(oldCor, digits=3), "->", round(newCor, digits=3)))

    addResult("ImputeAcc", filePrefix, round(newAcc, digits=3))
    addResult("ImputeYield", filePrefix, round(newYield, digits=3))
    addResult("ImputeCor", filePrefix, round(newCor, digits=3))

}
# assessImputation("genotypes.txt")
assessImputation("imputation.genotypes")


results = do.call("rbind", results)
write.table(results, output)



