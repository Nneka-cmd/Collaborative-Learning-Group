geno <- gl.load("original_recoced.rds")
karp <- gl.drop.pop(geno, pop.list = c('Rosenow', 'Quinby'), recalc = TRUE,mono.rm = TRUE)
quin <- gl.drop.pop(geno, pop.list = c('Rosenow', 'Karper'), recalc = TRUE,mono.rm = TRUE)
indNames(quin)
# pheno <- read.csv("karp_imputed.csv")
# pheno2 <- pheno[, c("id", "Origin")]
# pheno2$id <- as.factor(pheno2$id)
# pheno2$Classification <- as.factor(pheno2$Classification)
Fert_class <- Pheno2[, c("id", "class", "grain_vs_forage", "Genotype.1", "current_vs_historic")]

# Reorder Fert_class to match the tip labels
ordered2 <- Fert_class[match(indNames(quin), Fert_class$id),]

other(quin) <- as.list(ordered2)
pop(quin) <- quin$other$class
nInd(quin)
pcs <- gl.pcoa(quin)
gl.pcoa.plot(pcs, quin, pop.labels = "legend", ellipse = TRUE)

# DISTANCE MATRIX
D <- gl.dist.ind(karp)
pco <- gl.pcoa(D)
gl.pcoa.plot(pco,karp,ellipse=TRUE)

#write_csv(as.data.frame(indNames(karp)), "karp_filtered.txt")

# Extract eigenvalues and compute percent variance
eig_vals <- pcs$eig
pc_percent <- eig_vals / sum(eig_vals) * 100
# View first few PCs
round(pc_percent[1:5], 2)
# PC1: 25.31%, PC2: 13.47%, ...

library(ggplot2)

# Extract PCA scores from pcs
pcs_scores <- as.data.frame(pcs$scores)
pcs_scores$ID <- quin@ind.names  # Add sample IDs
pcs_scores$class <- quin@pop
pcs_scores$genotype <- quin@other$Genotype.1
pcs_scores$grain_forage <- quin@other$grain_vs_forage
pcs_scores$current_historic <- quin@other$current_vs_historic

# Plot with ggplot
ggplot(pcs_scores, aes(x = PC2, y = PC3, label = genotype, colour = current_historic)) +
  geom_point() +
  geom_text(size = 2, hjust = 0, vjust = 1.5) +
  theme_minimal()+
  labs(title = "Quinby Current vs Historic",
       x = paste0("PC1 (", round(pc_percent[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_percent[2], 1), "%)"))
