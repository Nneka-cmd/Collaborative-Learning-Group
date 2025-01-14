#Start working with individual RIL families from the NAM,
#rather the using the whole NAM population.
#You can assign QTL effects to marker directly in R for simple genetic
#architectures (e.g. see Bouchet et al. 2017) and/or use Lipka lab's
#simplePHENOTYPES for more sophisticated genetic architectures
# https://cran.r-project.org/web/packages/simplePHENOTYPES/index.html

rm(list = ls())
library(AlphaSimR)
library(lsa)

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

SP = SimParam$new(founders)
SP$setTrackPed(TRUE)

refLine = newPop(founders[1], id="R")
founderLines = newPop(founders[2:11], id=paste0("D", 1:10))

refGeno = pullSegSiteGeno(refLine)

crossPlan = matrix(c(rep(1:10), rep(1,10)), nrow=10, ncol=2)
F1 = makeCross2(founderLines, refLine, crossPlan, simParam=SP)
F2 = self(F1)
F3 = self(F2)
F4 = self(F3)
F5 = self(F4)
F6 = self(F5, nProgeny=nProgeny)

F6Geno = pullSegSiteGeno(F6)

# consider making explicit families of each pedigree
pca_founders <- prcomp(founderGeno)
eigen_founders <- pca_founders$sdev^2
genetic_distances_founders <- dist(founderGeno, method = "euclidean")
hclust_founders <- hclust(genetic_distances_founders, method = "average")

pca_F6 <- prcomp(F6Geno)
eigen_F6 <- pca_F6$sdev^2
genetic_distances_F6 <- dist(F6Geno, method = "euclidean")
hclust_F6 <- hclust(genetic_distances_F6, method = "average")

par(mfrow = c(3,2))
par(mar=c(1,1,1,1))
#color_palette <- rainbow(nFoundingInd) # numbers of colors to use
color_palette <- (c("blue", "red", "orange", "gold", "yellow", "darkorange3", "deeppink3", "khaki3", "lightsalmon", "tan3", "wheat2"))
colors <- rep(color_palette[2:11], each = nProgeny) # creating vector of colors

plot(eigen_founders, type = "b", col = "black",
     xlab = "Principal Component", ylab = "Eigenvalue", 
     main="Variance Explained PC in Founders")

plot(eigen_F6, type = "b", col = "black",
     xlab = "Principal Component", ylab = "Eigenvalue", 
     main="Variance Explained PC in RIL Families")

plot(x=pca_founders$x[,1],y=pca_founders$x[,2], main="PCA of Founders", col = color_palette)
plot(x=pca_F6$x[,1],y=pca_F6$x[,2], main="PCA of RIL Families", col = colors)

plot(hclust_founders, main = "Clustering of Founders", labels = FALSE)

#selected
plot(hclust_F6, main = "Clustering of RIL Families", labels = FALSE)
