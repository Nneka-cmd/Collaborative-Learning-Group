createRIL <- function(interPop=TRUE) {
  aIdx <- runif(2,1,nInd(popA))
  parentA <- popA[aIdx[1]]
  
  if (interPop) {
    parentB <- popB[runif(1,1,nInd(popB))]
  } else {
    parentB <- popA[aIdx[2]]
  }
  for (f in 1:10) {
    parentA <- self(parentA)
    parentB <- self(parentB)
  }
  
  F1 <- randCross2(parentA,
                   parentB,
                   nCrosses=1,
                   nProgeny=n.RILFams)
  
  F2 <- self(F1)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  RIL <- self(F9, nProgeny=n.indPerRILFam)
  #BC1 <- randCross2(parentA, F1, nCrosses=200)
  return (c(RIL, parentA, parentB))
}

createGWASPop <- function() {
  gwasPop <- c(popA, popB)
}

getABHGeno <- function(RIL, parentA, parentB, genoEnc=c(NA,1,2,3)) {
  # Filter out monomorphic snps
  allSnps <- as.data.frame(t(rbind(pullSnpGeno(parentA),pullSnpGeno(parentB))))
  colnames(allSnps) <- c('parentA', 'parentB')
  monoSnps <- allSnps %>% filter(parentA==parentB)
  monoMarkers <- rownames(monoSnps)
  
  aSnps <- pullSnpGeno(parentA) 
  bSnps <- pullSnpGeno(parentB) 
  rilSnps <- pullSnpGeno(RIL)
  output <- rilSnps
  genoEnc=c(NA,1,2,3)
  for (m in 1:ncol(rilSnps)) {
    if (colnames(rilSnps)[m] %in% monoMarkers) {
      output[,m] <- genoEnc[1]
      next
    }
    for (i in 1:nrow(rilSnps)) {
      rilGeno <- rilSnps[i,m]
      aGeno <- aSnps[1,m]
      bGeno <- bSnps[1,m]
      if ((aGeno == 1) || (bGeno == 1)) {
        output[i,m] <- genoEnc[1]
      } else if (rilGeno == 1) {
        output[i,m] <- genoEnc[3]
      } else if (rilGeno == aGeno) {
        output[i,m] <- genoEnc[2]
      } else if (rilGeno == bGeno) {
        output[i,m] <- genoEnc[4]
      }
    }
  }
  return (output)
}

getCross <- function(pop, parentA, parentB, popType) {
  snpGeno <- getABHGeno(pop, parentA, parentB)
  chrs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  pheno <- data.frame(pheno=pop@pheno)
  geno <- list()
  for (c in chrs) {
    snps <- snpGeno[,grepl(paste0(c,"_"), colnames(snpGeno))]
    chr_snps <- colnames(snps)
    chr_snps <- sapply(chr_snps, function(x) as.numeric(strsplit(x, split="_")[[1]][2])/10)
    chr_geno <- snps
    geno[[c]] <- list()
    attributes(geno[[c]])$class<-"A"
    geno[[c]]$"data" <- chr_geno
    geno[[c]]$"map" <- chr_snps
  }
  cross <- list(geno=geno,pheno=pheno)
  class(cross) <- c(popType,"cross")
  return (cross)
}