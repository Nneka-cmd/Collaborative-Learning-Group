#load packages
library(dartR.base)
library(tidyverse)
library(GenomicRanges)

#import data to R
All_companies <- gl.load("filtered2.rds")
#Confirming data has recoded individuals and populations

popNames(All_companies)
indNames(All_companies)

#Parse filtered data by companies
Quinby <- gl.drop.pop(All_companies, pop.list = c('Karper', 'Rosenow'),
                      mono.rm = TRUE, recalc = TRUE)
Karper <- gl.drop.pop(All_companies, pop.list = c('Quinby', 'Rosenow'),
                      mono.rm = TRUE, recalc = TRUE)
Rosenow <- gl.drop.pop(All_companies, pop.list = c('Quinby', 'Karper'),
                       mono.rm = TRUE, recalc = TRUE)

#Rosenow SNP density
rose_metrics <- as.data.frame(Rosenow@other$loc.metrics)
rose_subset <- rose_metrics %>%
  select(c("Chrom_Sorghum_bicolor_v5.1", "ChromPosSnp_Sorghum_bicolor_v5.1",
           "Strand_Sorghum_bicolor_v5.1")) %>%
  rename(chrom = "Chrom_Sorghum_bicolor_v5.1",
         position = "ChromPosSnp_Sorghum_bicolor_v5.1",
         strand = "Strand_Sorghum_bicolor_v5.1") %>%
  mutate(strand2 = case_when(strand ==  "Plus" ~ "+",
                             strand == "Minus" ~ "-"))
# Define bin size (e.g., 1 Mb) Windows
bin_size <- 1e6

# Compute SNP density by binning positions
snp_density <- snp_subset %>%
  group_by(chrom) %>%
  mutate(Bin = floor(position / bin_size) * bin_size) %>%
  count(chrom, Bin)

# Create SNP density plot
ggplot(rose_density, aes(x = Bin, y = n, fill = chrom)) +
  geom_bar(stat = "identity", width = bin_size * 0.8) +
  facet_wrap(~ chrom, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  labs(title = "SNP Density Across Genome",
       x = "Genomic Position (bp)",
       y = "SNP Count per Bin") +
  scale_fill_brewer(palette = "Set3")

# Define QTL positions (Replace with actual QTL data)
qtl_data <- read.csv("QTL positions - Sheet1.csv")

# Merge QTL data for plotting
qtl_plot_data <- qtl_data %>%
  mutate(y_min = 0, y_max = max(rose_density$n, na.rm = TRUE))  # Adjust QTL height

#Plot with snp densities with QTLs 
ggplot(snp_density, aes(x = Bin, y = n, fill = chrom)) +
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

#Confirming QTL postions on snp data
# Convert SNPs to GRanges object
snp_subset$chrom <- as.character(snp_subset$chrom)
snp_ranges <- GRanges(seqnames = snp_subset$chrom,
                       ranges = IRanges(start = snp_subset$position,
                                        end = snp_subset$position),
                      strand = snp_subset$strand2)
#convert genes to Granges
gene_ranges <- GRanges(seqnames = qtl_data$chrom,
                       ranges = IRanges(start = qtl_data$Start,
                                        end = qtl_data$End))
#Find overlaps (query snps against QLT)
overlaps <- findOverlaps(snp_ranges, gene_ranges)

#extract data overlap hits into a dataframe
gene_hits <- data_frame(
  chrom = as.character(seqnames(snp_ranges[queryHits(overlaps)])),
  snp_pos = start(snp_ranges[queryHits(overlaps)]))