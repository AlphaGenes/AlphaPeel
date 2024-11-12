rm(list = ls())
# install.packages(pkg = "AlphaSimR")
library(AlphaSimR)

# ---- Two (or more) sub-populations A and B which cross in later generations ----

parameters <- read.table("../simulation_parameters.txt")
nparams <- nrow(parameters)
for (parameter in (1:nparams)) {
  eval(parse(text = paste0(parameters$V1[parameter], "<-", parameters$V2[parameter])))
}

nIndPerGen <- nInd / nGen

nLociAllPerChr <- floor(nLociAll / nChr)
nLociHDPerChr <- nLociHD / nChr
nLociLDPerChr <- nLociLD / nChr

founderGenomes <- runMacs(nInd = nIndPerGen, nChr = nChr, segSites = nLociAllPerChr, split = 1000)
SP <- SimParam$new(founderPop = founderGenomes)
SP$setSexes("yes_rand")
SP$setTrackPed(isTrackPed = TRUE)
SP$setTrackRec(isTrackRec = TRUE)
varG <- matrix(data = c( 1.0, -0.3,
                         -0.3,  1.0), byrow = TRUE, nrow = 2)
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0), var = diag(varG), corA = varG)
varE <- matrix(data = c(2.0, 0.0,
                        0.0, 2.0), byrow = TRUE, nrow = 2)

collectData <- function(pop, data = NULL, population, generation) {
  remove <- FALSE
  if (is.null(data)) {
    remove <- TRUE
    data <- vector(mode = "list", length = 5)
    names(data) <- c("pedigree", "haploIBD", "genoIBS", "Haplotype", "Genotype")
    data$pedigree <- data.frame(id = NA, population = NA, generation = NA,
                                mid = NA, fid = NA,
                                gv1 = NA, pv1 = NA,
                                gv2 = NA, pv2 = NA)
    data$haploIBD <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$genoIBS <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$Haplotype <- matrix(data = NA, ncol = sum(pop@nLoci))
    data$Genotype <- matrix(data = NA, ncol = sum(pop@nLoci))
  }
  data$pedigree <- rbind(data$pedigree,
                         data.frame(id = pop@id,
                                    population = population,
                                    generation = generation,
                                    mid = pop@mother,
                                    fid = pop@father,
                                    gv1 = pop@gv[, 1],
                                    pv1 = pop@pheno[, 1],
                                    gv2 = pop@gv[, 2],
                                    pv2 = pop@pheno[, 2]))
  data$haploIBD <- rbind(data$haploIBD,
                         pullIbdHaplo(pop = pop))
  data$genoIBS <- rbind(data$genoIBS,
                        pullSegSiteGeno(pop = pop))
  listLoci <- seq.int(from = 1, to = sum(pop@nLoci), by = 1)
  marker <- ""
  i <- 1
  for (i in 1:length(listLoci)){
    ind <- paste("1_", listLoci[i], sep = "")
    marker <- c(marker, ind)
  }
  marker <- marker[-1]
  
  data$Haplotype <- rbind(data$Haplotype,
                         pullMarkerHaplo(pop = pop, markers = marker))
  data$Genotype <- rbind(data$Genotype,
                          pullMarkerGeno(pop = pop, markers = marker))
  if (remove) {
    data$pedigree <- data$pedigree[-1, ]
    data$haploIBD <- data$haploIBD[-1, ]
    data$genoIBS <- data$genoIBS[-1, ]
    data$Haplotype <- data$Haplotype[-1, ]
    data$Genotype <- data$Genotype[-1, ]
  }
  return(data)
}

# Founder population & split
founders <- newPop(rawPop = founderGenomes)
founders <- setPheno(pop = founders, varE = diag(varE))
popA <- founders[1:100]
popB <- founders[101:200]
data <- collectData(pop = popA, data = NULL, population = "A", generation = 0)
data <- collectData(pop = popB, data = data, population = "B", generation = 0)

# Select on each trait and keep the populations separate
for (generation in 1:nGen) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  if (generation == nGen){
    popA <- randCross(pop = parentsA, nCrosses = 75)
    popB <- randCross(pop = parentsB, nCrosses = 75)
  } else {
    popA <- randCross(pop = parentsA, nCrosses = 100)
    popB <- randCross(pop = parentsB, nCrosses = 100)
  }
  
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  data <- collectData(pop = popA, data = data, population = "A", generation = generation)
  data <- collectData(pop = popB, data = data, population = "B", generation = generation)
}

# Continued selection on each trait in each separate population,
# but add also continually admixed population selected on an index
popAB <- randCross(pop = c(parentsA, parentsB), nCrosses = 50)
popAB <- setPheno(pop = popAB, varE = diag(varE))
data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
economicWeights <- c(1, 1)
selIndexWeights <- smithHazel(econWt = economicWeights, varG = varG, varP = varG + varE)
for (generation in 6:10) {
  parentsA <- selectInd(pop = popA, nInd = 10, trait = 1)
  parentsB <- selectInd(pop = popB, nInd = 10, trait = 2)
  parentsAB <- selectInd(pop = popAB, nInd = 6, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsA4AB <- selectInd(pop = popA, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsB4AB <- selectInd(pop = popB, nInd = 2, trait = selIndex, scale = TRUE, b = selIndexWeights)
  parentsAB <- c(parentsAB, parentsA4AB, parentsB4AB)
  popA <- randCross(pop = parentsA, nCrosses = 75)
  popB <- randCross(pop = parentsB, nCrosses = 75)
  popAB <- randCross(pop = parentsAB, nCrosses = 50)
  popA <- setPheno(pop = popA, varE = diag(varE))
  popB <- setPheno(pop = popB, varE = diag(varE))
  popAB <- setPheno(pop = popAB, varE = diag(varE))
  data <- collectData(pop = popA,  data = data, population = "A",  generation = generation)
  data <- collectData(pop = popB,  data = data, population = "B",  generation = generation)
  data <- collectData(pop = popAB, data = data, population = "AB", generation = generation)
}

data$pedigree$population <- factor(data$pedigree$population, levels = c("A", "B", "AB"))
summary(data$pedigree$population)
data$pedigree$gv <- c(selIndex(Y = as.matrix(data$pedigree[, c("gv1", "gv2")]), scale = TRUE, b = selIndexWeights))
data$pedigree$pv <- c(selIndex(Y = as.matrix(data$pedigree[, c("pv1", "pv2")]), scale = TRUE, b = selIndexWeights))

data$pedigree$generationPlotShift <- data$pedigree$generation +
  c(-0.25, +0.25, 0)[data$pedigree$population]

# Get pedigree and assign metafounders
pedigree <- data.frame(data$pedigree$id, data$pedigree$mid, data$pedigree$fid, data$pedigree$population, data$pedigree$generation)
colnames(pedigree) <- c("id", "mid", "fid", "population", "generation")
pedigree <- pedigree[c(401:1400),]
pedigree$mid[c(1:100)] <- "MF_1"
pedigree$mid[c(101:200)] <- "MF_2"
pedigree$fid[c(1:100)] <- "MF_1"
pedigree$fid[c(101:200)] <- "MF_2"

pedigree <- pedigree[c(-4, -5)]
pedigree_noMF <- pedigree
pedigree_noMF$mid[c(1:200)] <- "0"
pedigree_noMF$fid[c(1:200)] <- "0"

# Get the genotypes (true) and observed (i.e with some missing)
genotypes <- data.frame(data$Genotype)
genotypes <- genotypes[c(401:1400),]

randFounder <- sort(sample.int(100, 70))
genotypes_obsv <- genotypes
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 70)) + 100
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 70)) + 600
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 70)) + 200
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 70)) + 700
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 40)) + 300
genotypes_obsv[randFounder, 1:2000] <- "9"
randFounder <- sort(sample.int(100, 40)) + 800
genotypes_obsv[randFounder, 1:2000] <- "9"

# Estimate the alternative allele frequencies.
founders_genotypes <- genotypes[c(1:100),]
alt_allele_prob_A <- colSums(founders_genotypes)
alt_allele_prob_A <- alt_allele_prob_A/(2*nrow(founders_genotypes))

founders_genotypes <- genotypes[c(101:200),]
alt_allele_prob_B <- colSums(founders_genotypes)
alt_allele_prob_B <- alt_allele_prob_B/(2*nrow(founders_genotypes))

alt_allele_prob_input_MF <- data.frame(matrix(nrow = 2, ncol = 2001))
alt_allele_prob_input_MF[1] <- c("MF_1", "MF_2")
alt_allele_prob_input_MF[1,2:2001] <- alt_allele_prob_A
alt_allele_prob_input_MF[2,2:2001] <- alt_allele_prob_B

founders_genotypes <- genotypes[c(1:200),]
alt_allele_prob <- colSums(founders_genotypes)
alt_allele_prob <- alt_allele_prob/(2*nrow(founders_genotypes))

alt_allele_prob_input_noMF <- data.frame(matrix(nrow = 1, ncol = 2001))
alt_allele_prob_input_noMF[1] <- "MF_1"
alt_allele_prob_input_noMF[1,2:2001] <- alt_allele_prob

# Haplotypes
# Get haplotypes for phased and unphased genotype probabilities
haplotypes <- data.frame(data$Haplotype)
haplotypes <- haplotypes[c(401:1400),]

for (ind in 1:nInd) {
  maternal <- haplotypes[ind * 2 - 1, ]
  haplotypes[ind * 2 - 1, ] <- haplotypes[ind * 2, ]
  haplotypes[ind * 2, ] <- maternal
}

# ----- Phased genotypes probability------
haplotypes <- data$Haplotype[c(801:2800),]

phasedGenotypes <- matrix(data = 0, nrow = nInd * 4, ncol = nLociAll + 1)
for (ind in (1:nInd)) {
  for (locus in (1:nLociAll)) {
    currentGeno <- haplotypes[((ind - 1) * 2 + 1):(ind * 2), locus]
    if (all(currentGeno == c(0, 0)) == TRUE) {
      currentPhasedGeno <- c(1, 0, 0, 0)
    } else if (all(currentGeno == c(0, 1)) == TRUE) {
      currentPhasedGeno<- c(0, 1, 0, 0)
    } else if (all(currentGeno == c(1, 0)) == TRUE) {
      currentPhasedGeno <- c(0, 0, 1, 0)
    } else {
      currentPhasedGeno <- c(0, 0, 0, 1)
    }
    phasedGenotypes[((ind - 1) * 4 + 1):(ind * 4), locus + 1] <- currentPhasedGeno
  }
  phasedGenotypes[((ind - 1) * 4 + 1):(ind * 4), 1] <- ind
}

# ----- (Unphased) genotype probability -----

UnphasedGenotypes <- matrix(data = 0, nrow = nInd * 3, ncol = nLociAll + 1)
for (ind in (1:nInd)) {
  for (locus in (1:nLociAll)) {
    currentGeno <- haplotypes[((ind - 1) * 2 + 1):(ind * 2), locus]
    if (all(currentGeno == c(0, 0)) == TRUE) {
      currentUnphasedGeno <- c(1, 0, 0)
    } else if (all(currentGeno == c(0, 1)) == TRUE) {
      currentUnphasedGeno<- c(0, 1, 0)
    } else if (all(currentGeno == c(1, 0)) == TRUE) {
      currentUnphasedGeno <- c(0, 1, 0)
    } else {
      currentUnphasedGeno <- c(0, 0, 1)
    }
    UnphasedGenotypes[((ind - 1) * 3 + 1):(ind * 3), locus + 1] <- currentUnphasedGeno
  }
  UnphasedGenotypes[((ind - 1) * 3 + 1):(ind * 3), 1] <- ind
}

# To match the ids (+ 400)
phasedGenotypes[,1] <- phasedGenotypes[,1] + 400
UnphasedGenotypes[,1] <- UnphasedGenotypes[,1] + 400

# Save the results/data for testing
write.table(x = pedigree, file = "metafounder_ped_file.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = genotypes_obsv, file = "metafounder_geno_file.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(x = genotypes, file = "true-metafounder_dosage.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(x = UnphasedGenotypes, file = "true-metafounder_geno_prob.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = phasedGenotypes, file = "true-metafounder_phased_geno_prob.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = alt_allele_prob_input_MF, file = "metafounder_alt_allele_prob_file.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


