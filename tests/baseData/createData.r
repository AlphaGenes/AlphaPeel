

trueGenotypes = read.table("trueGenotypes.txt")
pedigree = read.table("pedigree.txt")

#Alternative pedigrees for assign.

pedigreeNoParents = pedigree
pedigreeNoParents[,c(2,3)] = 0
write.table(pedigreeNoParents, "pedigree.noParents", row.names=F, col.names=F, quote=F)

pedigreeNoDam = pedigree
pedigreeNoDam[,3] = 0
write.table(pedigreeNoParents, "pedigree.noDam", row.names=F, col.names=F, quote=F)

lastGen = 801:1000
sireList = unique(pedigree[lastGen,2])

sireList = do.call("rbind", lapply(lastGen, function(id) {c(id, sireList)}))

write.table(sireList, "sire.list", row.names=F, col.names=F, quote=F)

#Create stuff for MGS Assign

grandInfo = lapply(1:nrow(pedigree), function(index) {
        ind = pedigree[index,1]
        mother = pedigree[index,3]
        if(mother == 0) return(c(0,0))

        motherIndex = which(pedigree[,1] == mother)
        return(pedigree[motherIndex, c(2,3)])
    })
grandInfo = do.call("rbind", grandInfo)

pedigreeWithMatGrand = data.frame(pedigree[,1], pedigree[,2], pedigree[,3], grandInfo[,1], grandInfo[,2], stringsAsFactors = FALSE)
# pedigree[,3] = 0


grandSires = unique(pedigreeWithMatGrand[lastGen, 4])
grandSireList = t(sapply(lastGen, function(id) {return(c(pedigree[id,1], grandSires))}))

write.table(grandSireList, "grandsire.list", row.names=F, col.names=F, quote=F)


write.table(pedigreeWithMatGrand, "pedigree.extended", row.names=F, col.names=F, quote=F)

### Map files for peeling

nLoci = ncol(trueGenotypes) -1
values = data.frame(1, paste0("1-", 1:nLoci), 1:nLoci) 
write.table(values, "map.txt", row.names=F, col.names=F, quote=F)

subset = floor(seq(1, nLoci, length.out = 200))
subsetValues = values[subset,]
write.table(subsetValues, "segmap.txt", row.names=F, col.names=F, quote=F)


### Sequence data

sequence = c()
nLoci = ncol(trueGenotypes) -1
coverage = 2
gbsRandom = function(pos) {

    id = trueGenotypes[pos,1]
    geno = trueGenotypes[pos, ]

    selectedLoci = 1:nLoci

    # loci = sample(selectedLoci, coverage, replace=T)
    # counts = rep(0, nLoci)
    # for( val in loci){
    #     counts[val] = counts[val] + 1
    # }
    counts = rpois(nLoci, coverage)
    error = 0.001
    probs = (1-error) * as.numeric(geno[,-1])/2.0 + error/2
    
    probs = probs*sign(counts) #This will fix any missing values IF the loci is not selected. Will not fix and will produce NA if the loci is selected.

    altReads = rbinom(nLoci, counts, prob = probs)

    values = matrix(c(id, counts-altReads, id, altReads), nrow = 2, byrow=TRUE)

    if(length(sequence) == 0) {
        sequence <<- values
    }
    else{
        sequence <<- rbind(sequence, values)
    }

}

for(i in 1:nrow(trueGenotypes)) {
    gbsRandom(i)
}
write.table(sequence, "sequence.2", row.names=F, col.names=F, quote=F)


