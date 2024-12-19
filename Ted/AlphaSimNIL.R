# Now that you've made mock data sets and visual predictions for a simple NIL experiment, a next step would be to add realistic complexity to the experimental design and predictions.
# Define a specific complex trait (e.g. biotic resistance, abiotic tolerance, or nutritional value)
# that you're targeting and an experimental approach that would be required to phenotype that trait
# What controls would be needed to demonstrate:
#   That the phenotype is what you're purporting it to be, and not some other experimental artifact
#   That the level of the trait difference (allelic substitution effect) is large enough to be
#     relevant to breeding
#   The NILs have expected phenotypes relative to the recurrent and donor parent
#   You'll want to have both positive and negative controls
# GxE interactions: https://link.springer.com/chapter/10.1007/978-981-13-7095-3_20

# THIS IS A WORK-IN-PROGRESS

library(AlphaSimR)
library(lsa)
rm(list = ls())
par(mfrow = c(1, 1))


nQTLs = 100
nChromosomes = 11
nProgeny = 10
nFoundingInd = 11

founders = runMacs(
  nInd=nFoundingInd,
  nChr=nChromosomes,
  segSites=nQTLs,
  species = "MAIZE",
  inbred = TRUE,
)
founderGeno = pullSegSiteGeno(founders)

nilSP$addTraitA(
  mean=4.5,
  var=1,
  nQtlPerChr=5,
)$setVarE(H2 = 0.5)



# Use importInbredGeno for creating NILs
#https://jvanderw.une.edu.au/03Practical_Breeding%20program/Extra/traitIntrogression.R
nMarkers = 11
nChr = 2

genMap = data.frame(
  markerName = paste0("M", 1:(nChr * nMarkers)),
  chromosome = rep(1:nChr, each = nMarkers),
  position = rep(seq(from = 0, to = 1, by = 0.1), times = 2)
)

geno = matrix(
  sample(2:0, nMarkers*nChr, TRUE),
  nrow = nChr, ncol = nMarkers
)

colnames(geno) = genMap$markerName

ped = data.frame(
  id = c("RP", "Donor"),
  mother = c(0, 0),
  father = c(0, 0)
)


founderPop = importInbredGeno(geno = geno,
                              genMap = genMap,
                              ped = ped)
nilSP = SimParam$new(nilFounders)

nilSP$addTraitA(
  mean=120,
  var=100,
  nQtlPerChr=5,
)$setVarE(H2 = 0.5)

nilSP$traits