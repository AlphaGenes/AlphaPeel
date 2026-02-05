args = commandArgs(trailingOnly = TRUE)
test = args[1]


readAndSortFile = function(fileName, idList = NULL, ...) {

    values = read.table(fileName, ..., stringsAsFactors = FALSE)
    if(!is.null(idList)) {
        values = values[values[,1] %in% idList,]
    }
    values = values[order(values[,1]),]
    return(values)

}

if(test == "1") {

    true = readAndSortFile("test1/trueGenotypes.txt")
    values = readAndSortFile("test1/outputs/output.called.0.1")

    if(all(true == values)) {
        print("Test 1: Success")
    }else{
        print("Test 1: TEST 1 FAILED")

    }
}



if(test == "2") {

    true = readAndSortFile("test2/trueGenotypes.txt")[,c(1, 3, 4, 5)]
    values = readAndSortFile("test2/outputs/output.called.0.1")

    if(all(true == values)) {
        print("Test 2: Success")
    }else{
        print("Test 2: TEST 2 FAILED")
        print(true)
        print(values)

    }
}


if(test == "3") {

    true = readAndSortFile("test3/trueGenotypes.txt")[,c(1, 3, 4, 5)]
    #This is in plink format
    values = readAndSortFile("test3/outputs/output.called.0.99.raw", header= TRUE)[, c(2, 7, 8, 9)]
    values = values[order(values[,1]),]
    values[is.na(values)] = 9

    success = TRUE
    if(!all(true == values)) {
        success = FALSE
        print("Test 3: Failed on .raw file")
    }

    bimFile = read.table("test3/outputs/output.called.0.99.bim")
    target = data.frame(c(1, 1, 1), c("snp2", "snp3", "snp4"), c(0, 0, 0), c(2, 3, 4), c("A", "A", "A"), c("B", "B", "B"))

    if(!all(bimFile == target)){
        success = FALSE
        print("Test 3: Failed on .bim file")
    }

    if(success) {
        print("Test3: raw file success")
        print("Test3: bim file success")
        print("Test 3: Success")
    }else{
        print("Test 3: TEST 3 FAILED")
        print(true)
        print(values)

    }
}


if(test == "3b") {
    success = TRUE

    for(nind in 1:4) {
        true = readAndSortFile(paste0("test3b/outputs/output.", nind, ".called.0.1.raw"), header= TRUE)
        values = readAndSortFile(paste0("test3b/outputs/output.", nind, ".called.0.1.raw"), header= TRUE)
        if(!all(true == values)) {
            success = FALSE
            print(paste("Test 3: Failed on nInd", nind))
        }   

    }

    if(success) {
        print("Test 3b: Success")
    }else{
        print("Test 3b: TEST 3 FAILED")
    }
}



if(test == "3c") {
    success = TRUE

    true = readAndSortFile("test3c/genotypes.txt")
    values = readAndSortFile("test3c/outputs/noFamNoPedigree.called.0.9")
    if(!all(true == values)) {
        success = FALSE
        print("Test 3c: Failed on noFamNoPedigree")
    }   

    true = readAndSortFile("test3c/trueGenotypes.txt")
    values = readAndSortFile("test3c/outputs/noFamPedigree.called.0.9")
    if(!all(true == values)) {
        success = FALSE
        print("Test 3c: Failed on noFamPedigree.called.0.9")
    }   

    true = readAndSortFile("test3c/trueGenotypes.txt")
    values = readAndSortFile("test3c/outputs/famNoPedigree.called.0.9")
    if(!all(true == values)) {
        success = FALSE
        print("Test 3c: Failed on famNoPedigree.called.0.9")
    }   


    if(success) {
        print("Test 3c: Success")
    }else{
        print("Test 3c: TEST 3 FAILED")
    }
}


if(test == "4") {
    methods = c("id","pedigree","genotypes","sequence")
    answer = c("genotypes","penetrance","genotypes","seq")
    for(i in 1:length(methods)) {
        values = read.table(paste0("test4/outputs/output.", methods[i], ".called.0.1"))

        if(nrow(values) == 4 & values[1, 1] == answer[i]) {
            print(paste0("Test 4a ", methods[i], ": Success"))
        }else{
            print(paste0("Test 4a ", methods[i], ": FAILED"))
            print(values)
        }
    }
    values = read.table("test4/outputs/output.only.called.0.1")
    if(nrow(values) == 1 & values[1, 1] == "seq") {
        print("Test 4b: Success")
    }else{
        print("Test 4b: TEST 2 FAILED")
        print(values)
    }
}


checkForFiles = function(filePrefix) {

    result = c( file.exists(paste0(filePrefix, ".dosages")),
                file.exists(paste0(filePrefix, ".seg")),
                file.exists(paste0(filePrefix, ".maf")),
                file.exists(paste0(filePrefix, ".genoError")),
                file.exists(paste0(filePrefix, ".seqError")),
                file.exists(paste0(filePrefix, ".haps"))
     )
    return(result)
}



if(test == "5") {
    print("TEST 5 NOT IMPLIMENTED!!!!!")
}



if(test == "7") {

    true = readAndSortFile("test7/trueGenotypes.txt")
    values = readAndSortFile("test7/outputs/output.called.0.1")

    if(all(true == values)) {
        print("Test 7: Success")
    }else{
        print("Test 7: TEST 7 FAILED")

    }
}


if(test == "7b") {

    true = readAndSortFile("test7b/trueGenotypes.txt")
    values = readAndSortFile("test7b/outputs/output.called.0.1")

    if(all(true == values)) {
        print("Test 7b: Success")
    }else{
        print("Test 7b: TEST 7 FAILED")

    }
}


if(test == "7c") {

    true = readAndSortFile("test7c/trueGenotypes.txt")
    values = readAndSortFile("test7c/outputs/output.called.0.1")

    if(all(true == values)) {
        print("Test 7c: Success")
    }else{
        print("Test 7c: TEST 7 FAILED")

    }
}



if(test == "8") {
    no_dosage = all(checkForFiles("test8/outputs/no_dosage") == c(F,T, T, T, T, F))
    no_seg = all(checkForFiles("test8/outputs/no_seg") == c(T, F, T, T, T, F))
    no_params = all(checkForFiles("test8/outputs/no_params") == c(T, T, F, F, F, F))
    haps = all(checkForFiles("test8/outputs/haps") == c(T, T, T, T, T, T))

    success = TRUE
    if(!no_dosage){
        print("Test 8: no_dosage failed")
        success = FALSE
    }
    if(!no_seg){
        print("Test 8: no_seg failed")
        success = FALSE
    }
    if(!no_params){
        print("Test 8: no_params failed")
        success = FALSE
    }
    if(!haps){
        print("Test 8: haps failed")
        success = FALSE
    }

    if(success) {
        print("Test 8: Success!")
    }

}


