library(tidyverse)
library(ggrepel)
rm(dist_matrix, pop_group, gi)

all_data <- gl.load("all_imputed.rds")
nLoc(all_data)
locNames(all_data)

rose <- gl.drop.pop(all_data, pop.list = c('Karp', 'Quin'), recalc = TRUE,mono.rm = TRUE)
quin <- gl.drop.pop(all_data, pop.list = c('Karp', 'Rose'), recalc = TRUE,mono.rm = TRUE)
karp <- gl.drop.pop(all_data, pop.list = c('Rose', 'Quin'), recalc = TRUE,mono.rm = TRUE)

rose_subset <- as.data.frame(cbind(rose@chromosome, rose@position)) %>%
  rename(chrom = "V1", position = "V2")
rose_subset$chrom <- as.factor(rose_subset$chrom)

quin_subset <- as.data.frame(cbind(quin@chromosome, quin@position)) %>%
  rename(chrom = "V1", position = "V2")
quin_subset$chrom <- as.factor(quin_subset$chrom)

karp_subset <- as.data.frame(cbind(karp@chromosome, karp@position)) %>%
  rename(chrom = "V1", position = "V2")
karp_subset$chrom <- as.factor(karp_subset$chrom)

# Define bin size (e.g., 1 Mb) Windows
bin_size <- 1e6

qtl_data <- read.csv("QTL positions - Sheet1.csv")
qtl_data$chrom <- as.factor(qtl_data$chrom)

# Compute SNP density by binning positions
rose_density <- rose_subset %>%
  group_by(chrom) %>%
  mutate(Bin = floor(position / bin_size) * bin_size) %>%
  count(chrom, Bin)
# Merge QTL data for plotting
qtl_plot_data <- qtl_data %>%
  mutate(y_min = 0, y_max = max(rose_density$n, na.rm = TRUE))
#Plot with snp densities with QTLs 
ggplot(rose_density, aes(x = Bin, y = n, fill = chrom)) +
  geom_bar(stat = "identity", width = bin_size * 0.8) +
  geom_rect(data = qtl_plot_data,
            aes(xmin = Gene_Start,xmax = Gene_End,
                ymin = y_min, ymax = y_max, fill = chrom),
            color = "black", alpha = 0.3, inherit.aes = FALSE) +
  geom_text_repel(data = qtl_plot_data,
                  aes(x = (Gene_Start + Gene_End) / 2,  # Midpoint of QTL region
                      y = y_max * 1.02,  # Slightly below the max SNP count
                      label = Gene.Nomin),
                  hjust = 0,color = "black", size = 4,
                  fontface = "bold", inherit.aes = FALSE,
                  nudge_x = 10,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 2)+
  facet_wrap(~ chrom, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Rosenow SNP Density with classical genes",
       x = "Genomic Position (bp)",
       y = "SNP Count per Bin")

# Compute SNP density by binning positions
quin_density <- quin_subset %>%
  group_by(chrom) %>%
  mutate(Bin = floor(position / bin_size) * bin_size) %>%
  count(chrom, Bin)
# Merge QTL data for plotting
qtl_plot_data <- qtl_data %>%
  mutate(y_min = 0, y_max = max(quin_density$n, na.rm = TRUE))
#Plot with snp densities with QTLs 
ggplot(quin_density, aes(x = Bin, y = n, fill = chrom)) +
  geom_bar(stat = "identity", width = bin_size * 0.8) +
  geom_rect(data = qtl_plot_data,
            aes(xmin = Gene_Start,xmax = Gene_End,
                ymin = y_min, ymax = y_max, fill = chrom),
            color = "black", alpha = 0.3, inherit.aes = FALSE) +
  geom_text_repel(data = qtl_plot_data,
                  aes(x = (Gene_Start + Gene_End) / 2,  # Midpoint of QTL region
                      y = y_max * 1.02,  # Slightly below the max SNP count
                      label = Gene.Nomin),
                  hjust = 0,color = "black", size = 4,
                  fontface = "bold", inherit.aes = FALSE,
                  nudge_x = 10,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 2)+
  facet_wrap(~ chrom, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Quinby SNP Density with classical genes",
       x = "Genomic Position (bp)",
       y = "SNP Count per Bin")
# Compute SNP density by binning positions
karp_density <- karp_subset %>%
  group_by(chrom) %>%
  mutate(Bin = floor(position / bin_size) * bin_size) %>%
  count(chrom, Bin)
# Merge QTL data for plotting
qtl_plot_data <- qtl_data %>%
  mutate(y_min = 0, y_max = max(karp_density$n, na.rm = TRUE))
#Plot with snp densities with QTLs 
ggplot(karp_density, aes(x = Bin, y = n, fill = chrom)) +
  geom_bar(stat = "identity", width = bin_size * 0.8) +
  geom_rect(data = qtl_plot_data,
            aes(xmin = Gene_Start,xmax = Gene_End,
                ymin = y_min, ymax = y_max, fill = chrom),
            color = "black", alpha = 0.3, inherit.aes = FALSE) +
  geom_text_repel(data = qtl_plot_data,
                  aes(x = (Gene_Start + Gene_End) / 2,  # Midpoint of QTL region
                      y = y_max * 1.02,  # Slightly below the max SNP count
                      label = Gene.Nomin),
                  hjust = 0,color = "black", size = 4,
                  fontface = "bold", inherit.aes = FALSE,
                  nudge_x = 10,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 2)+
  facet_wrap(~ chrom, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Karper SNP Density with classical genes",
       x = "Genomic Position (bp)",
       y = "SNP Count per Bin")




