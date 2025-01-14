# Fitness Functions

nChr <- 10

oneTraitFitFunc <- function(x) {
  return (-(x)^2) + rnorm(1,sd=2)
}

twoTraitFitFunc <- function(x) {
  res <- -((x[,1])^2) - ((x[,2])^2)
  return (res)
}

calculateFitnessTwoTrait <- function(x,y) {
  res <- -((x)^2) -(y^2)
  return (res)
}

hetLocus <- function(locus) {
  return (length(unique(locus)) > 1)
}

getEffectSize <- function(locus,
                          id,
                          fitFunc,
                          pop,
                          methodType) {
  
  if (!hetLocus(locus)) {
    return (NA)
  }
  
  strs = unlist(strsplit(id, "_"))
  chr = strtoi(strs[1])
  site = strtoi(strs[2])
  popSize <- nInd(pop)
  fitPre <- mean(fitFunc(gv(pop)))
  pop1 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=0, simParam=SP)
  pop2 <- editGenome(pop, ind=c(1:popSize), chr=chr, segSites=site, allele=1, simParam=SP)
  if (methodType == "MethodA") {
    effect_size_1 <- abs(mean(fitFunc(gv(pop1))) - fitPre)
    effect_size_2 <- abs(mean(fitFunc(gv(pop2))) - fitPre)
    return (max(effect_size_1, effect_size_2))
  } else if (methodType == "MethodB") {
    return (abs(mean(fitFunc(gv(pop2))) - mean(fitFunc(gv(pop1)))))
  }
}