set.seed(123)

basicPop = TRUE

if (basicPop) {
  founders = runMacs(
    nInd=n.popSize,
    nChr=n.chr,
    segSites=n.segSites,
    nThreads=4
  )
} else {
  founders = runMacs(nInd=n.ne,
                     nChr=10,
                     segSites=n.segSites,
                     inbred=TRUE,
                     manualCommand = paste(
                       "1000000000 -t", #Physical length 1e8 base pairs
                       2.5/1E8*(4*n.ne), #Mutation rate adjusted for Ne
                       "-r",1/1E8*(4*n.ne), #Recombination rate adjusted for Ne
                       "-eN",10/(4*n.ne),100/n.ne), #Modeling Ne=100 at 10 generations ago
                     manualGenLen=1, nThreads=4) #Genetic length 1 Morgan
}


SP <- SimParam$new(founders)
SP$restrSegSites(overlap = T)
#SP$addTraitA(mean=c(runif(1,-10,10),runif(1,-10,10)), var=c(1,1), nQtlPerChr=20)

SP$addTraitA(mean=n.initTraitVal, var=n.rate, nQtlPerChr=n.qtlPerChr, gamma = TRUE, shape = n.shape) # e.g. height
SP$addTraitA(mean=n.initTraitVal, var=n.rate, nQtlPerChr=n.qtlPerChr, gamma = TRUE, shape = n.shape) # e.g. flowering time
#SP$addTraitA(mean=950, var=200, nQtlPerChr = 100, gamma = TRUE, shape = 1)

SP$setVarE(h2=c(n.h2, n.h2))
SP$addSnpChip(nSnpPerChr=n.markers)
foundingPop <- newPop(founders, simParam = SP)

#plot3dPopulationFitness(foundingPop, calculateFitnessTwoTrait)
#plotTraitArchitecture(pop=foundingPop, method="MethodB", fitFunc=twoTraitFitFunc())