rm(list = ls())
library(AlphaSimR)
library(Rlab)

nLoci <- 1000
nInd <- 1000
nGen <- 5
nChr <- 1
nsegSites <- nLoci/nChr
genoError <- 0.001
seqError <- 0.01
nSelect <- 100

# runMacs for generic species
founderGenomes <- runMacs(nInd = (nInd / nGen),
                          nChr = nChr,
                          segSites = nsegSites,
                          species = "GENERIC")
SP <- SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
basePop <- newPop(founderGenomes)

currentPop <- basePop
Parents <- c()
for (generation in 1:(nGen - 1)) {
  # add selections, pick 100 parents randomly, make the parents high density, everyone in gen5 low density
  newParents <- selectInd(currentPop, nInd = nSelect, use = "rand")
  Parents <- c(Parents, getPed(newParents)[, 1])
  newPop <- randCross(currentPop, nCrosses = nInd(currentPop))
  currentPop <- newPop
  popList <- list(basePop, currentPop)
  basePop <- mergePops(popList)
}

# for the calculation of segregation and the recombination rate
# recHist[[individual]][[chromosome]][[haplotype]] 
# column 1 gives the parental haplotype, column 2 gives the start loci
recHist <- SP$recHist
IbdHaplo <- pullIbdHaplo(basePop) # flip maternal and paternal haplotype before use

# generate errors in the genotype input
# locus specific error rate?
generateGenoErr <- function(geno, error = genoError) {
  for (ind in 1:nInd){
    error_ind <- rbern(nLoci, error)
    for (locus in 1:nLoci){
      if (error_ind[locus] == 1) {
        if (geno[ind, locus] == 0) {
          geno[ind, locus] <- sample(c(1, 2), size = 1)
        }
        else {
          if (geno[ind, locus] == 1) {
            geno[ind, locus] <- sample(c(0, 2), size = 1)
          }
          else {
            geno[ind, locus] <- sample(c(0, 1), size = 1)
          }
        }
      }
    }
  }
  return(geno)
}

# generating genotypes
genotypes <- pullSegSiteGeno(basePop)

# generating haplotypes
haplotypes <- pullSegSiteHaplo(basePop)
# flip maternal and the paternal haplotype
for (ind in 1:nInd) {
  maternal <- haplotypes[ind * 2 - 1, ]
  haplotypes[ind * 2 - 1, ] <- haplotypes[ind * 2, ]
  haplotypes[ind * 2, ] <- maternal
}

url <- "https://gist.githubusercontent.com/gregorgorjanc/e7ca7a02ed59242573c12a890cf8c871/raw/e88ad2a3036f3385a0724d9d92a7cf5012ee0f86/simulateSeqReads.R"
import::from(pins::pin(url), "simulateSeqReads", .character_only=TRUE)

# generating sequence reads
sequence <- simulateSeqReads(genotypes, depth = 3, error = seqError)

# generating pedigree
pedigree <- getPed(basePop)

# generating genotypes with errors
genotypes_w_error <- generateGenoErr(genotypes, error = genoError)

# generating genotypes with missing data
# systematic missing 
# low density with 100 markers (1 marker known and 9 unknowns)
geno_missing <- matrix(rep(c(rep(1, times = 9), 0), times = nLoci / 10 * nInd), nrow = nInd, ncol = nLoci, byrow = TRUE)
for (ind in Parents) {
  geno_missing[as.numeric(ind), ] <- matrix(0, nrow = 1, ncol = nLoci)
}
genotypes_w_missing <- genotypes_w_error
for (ind in 1:nInd) {
  genotypes_w_missing[ind, geno_missing[ind, ] == 1] <- 9
}

# generating genotypes error rate
geno_error <- matrix(0, nrow = nLoci, ncol = 1)
for (locus in 1:nLoci) {
  total <- 0
  nerrors <- 0
  for (ind in 1:nInd) {
    if (genotypes_w_missing[ind, locus] != 9) {
      total <- total + 1
      if (genotypes_w_error[ind, locus] != genotypes[ind, locus]) {
        nerrors <- nerrors + 1
      }
    }
  }
  geno_error[locus, ] <- nerrors / total
}

# generating sequence error rate
seq_error <- matrix(0, nrow = nLoci, ncol = 1)
n_ref <- matrix(0, nrow = nInd, ncol = 1)
n_alt <- matrix(0, nrow = nInd, ncol = 1)

for (locus in 1:nLoci) {
  p_aa <- matrix(0, nrow = nInd, ncol = 1)
  p_AA <- matrix(0, nrow = nInd, ncol = 1)
  genotype_locus = genotypes[, locus]
  for (ind in 1:nInd) {
    if (genotype_locus[ind] == 0) {
      p_aa[ind, ] <- 1
    }
    else {
      if (genotype_locus[ind] == 2) {
        p_AA[ind, ] <- 1
      }
    }
    n_ref[ind, ] <- sequence[(ind - 1) * 2 + 1, locus]
    n_alt[ind, ] <- sequence[ind * 2, locus]
  }
  seq_error[locus, ] <- sum(n_alt * p_aa + n_ref * p_AA) / sum((n_ref + n_alt) * (p_aa + p_AA))
}

# generating minor allele frequency
maf <- matrix(0, nrow = nInd, ncol = 1)

for (locus in 1:nLoci) {
  count <- 0
  for (ind in 1:nInd) {
    count <- count + genotypes[ind, locus]
  }
  maf[locus, ] <- count / (nInd * 2)
}

# generating segregation file
segregation <- matrix(0, nrow = nInd * 4, ncol = nLoci + 1)
segregation[1:(4 * nInd / nGen), 2:(nLoci + 1)] <- matrix(rep(0.25, times = nLoci * 4 * nInd / nGen), nrow = 4 * nInd / nGen, ncol = nLoci)

for (ind in ((nInd / nGen) + 1):nInd) {
  indrecHist <- recHist[[ind]][[1]]
  
  nMaternalComb <- nrow(indrecHist[[1]])
  maternalPattern <- matrix(0, nrow = 4, ncol = nLoci)
  for (comb in (1:nMaternalComb)) {
    start <- indrecHist[[1]][comb, 2]
    if (indrecHist[[1]][comb, 1] == 1) {
      # maternal haplotype is inherited from the individual's grandmaternal
      pattern <- c(0, 1, 0, 1)
    }
    else {
      # maternal haplotype is inherited from the individual's grandpaternal
      pattern <- c(1, 0, 1, 0)
    }
    maternalPattern[1:4, start:nLoci] <- matrix(rep(pattern, times = nLoci - start + 1), nrow = 4, ncol = nLoci - start + 1)
  }
  
  nPaternalComb <- nrow(indrecHist[[2]])
  paternalPattern <- matrix(0, nrow = 4, ncol = nLoci)
  for (comb in (1:nPaternalComb)) {
    start <- indrecHist[[2]][comb, 2]
    if (indrecHist[[2]][comb, 1] == 1) {
      # paternal haplotype is inherited from the individual's grandmaternal
      pattern <- c(0, 0, 1, 1)
    }
    else {
      # paternal haplotype is inherited from the individual's grandpaternal
      pattern <- c(1, 1, 0, 0)
    }
    paternalPattern[1:4, start:nLoci] <- matrix(rep(pattern, times = nLoci - start + 1), nrow = 4, ncol = nLoci - start + 1)
  }
  
  startRow <- 4 * (ind - 1) + 1
  endRow <- 4 * ind
  segregation[startRow:endRow, 1] <- ind
  segregation[startRow:endRow, 2 : (nLoci + 1)] <- maternalPattern * paternalPattern
}

# generating genotype probability file
geno_prob <- matrix(0, nrow = nInd * 4, ncol = nLoci + 1)
for (ind in (1:nInd)) {
  for (locus in (1:nLoci)) {
    currentGeno <- haplotypes[((ind - 1) * 2 + 1):(ind * 2), locus]
    if (all(currentGeno == c(0, 0)) == TRUE) {
      current_geno_prob <- c(1, 0, 0, 0)
    }
    else {
      if (all(currentGeno == c(0, 1)) == TRUE) {
        current_geno_prob <- c(0, 1, 0, 0)
      }
      else {
        if (all(currentGeno == c(1, 0)) == TRUE) {
          current_geno_prob <- c(0, 0, 1, 0)
        }
        else {
          current_geno_prob <- c(0, 0, 0, 1)
        }
      }
    }
    geno_prob[((ind - 1) * 4 + 1):(ind * 4), 1] <- ind
    geno_prob[((ind - 1) * 4 + 1):(ind * 4), (locus + 1)] <- matrix(current_geno_prob, nrow = 4, ncol = 1)
  }
}

# TODO: generating recombination rate

# write the files
# modify the name for genotypes and haplotypes files later when calling threshold is known
write.table(genotypes, "true-called.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(haplotypes, "true-called_phase.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(sequence, "seqfile.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pedigree, "pedigree.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(geno_error, "true-genoError.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(seq_error, "true-seqError.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(maf, "true-maf.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genotypes_w_missing, "genotypes.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(segregation, "true-seg.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(geno_prob, "true-haps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Map files for peeling

nLoci = ncol(genotypes)
values = data.frame(1, paste0("1-", 1:nLoci), 1:nLoci)
write.table(values, "map.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

subset = floor(seq(1, nLoci, length.out = 200))
subsetValues = values[subset,]
write.table(subsetValues, "segmap.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

