createIndependentPops <- function() {
  pop <- foundingPop
  
  # BURN-IN
  for (gen in 1:n.burnInGens) {
    pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.burnInSelProp, nCrosses=nInd(pop))
  }
  
  # Create 2 sub populations
  # Figure out a better way to select randomly? or by traitA, traitB?
  popA <- selectInd(pop, use="rand", nInd=n.subPopSize)
  popB <- selectInd(pop, use="rand", nInd=n.subPopSize)
  
  for (gen in 1:n.gens) {
    if (mean(twoTraitFitFunc(pheno(popA))) < n.margin) {
      popA <- selectCross(popA, trait=twoTraitFitFunc, nInd=nInd(popA)*n.selProp, nCrosses=nInd(popA))
    }
    if (mean(twoTraitFitFunc(pheno(popB))) < n.margin) {
      popB <- selectCross(popB, trait=twoTraitFitFunc, nInd=nInd(popB)*n.selProp, nCrosses=nInd(popB))
    }
  }
  return (c(popA, popB))
}