#install R-CMplot for SM+NP density plot
install.packages("CMplot")
library("CMplot")

# if you want to use the latest version on GitHub:
source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")

Chromosome Haplotype Block Count
# Load necessary library
library(ggplot2)

# Create a data frame with genetic distance and R² values for each population
genetic_distance <- seq(0, 1, length.out = 100)
ld_magic <- exp(-5 * genetic_distance)       # MAGIC: Fast decay
ld_assoc <- exp(-3 * genetic_distance)       # Association panel: Moderate decay
ld_biparental <- exp(-1 * genetic_distance)  # Biparental: Slow decay

# Combine data into a data frame
ld_data <- data.frame(
  Genetic_Distance = rep(genetic_distance, 3),
  LD = c(ld_magic, ld_assoc, ld_biparental),
  Population = factor(rep(c("Biparental", "Association Panel", "GORP"), each = 100))
)

# Create the plot
ggplot(ld_data, aes(x = Genetic_Distance, y = LD, color = Population)) +
  geom_line(size = 1.2) +
  labs(
    title = "LD Decay (R²) vs Genetic Distance",
    x = "Genetic Distance (cM)",
    y = "Linkage Disequilibrium (R²)"
  ) +
  scale_color_manual(values = c("blue", "green", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_blank()
  )
library(tidyverse)

# Load necessary library
library(ggplot2)

# Create a data frame with genetic distance and R² values for each population
genetic_distance <- seq(0, 1, length.out = 100)
ld_assoc <- exp(-5 * genetic_distance)       # Association panel: Fast decay
ld_magic <- exp(-3 * genetic_distance)       # MAGIC: Moderate decay
ld_biparental <- exp(-1 * genetic_distance)  # Biparental: Slow decay

# Combine data into a data frame
ld_data <- data.frame(
  Genetic_Distance = rep(genetic_distance, 3),
  LD = c(ld_assoc, ld_magic, ld_biparental),
  Population = factor(rep(c("Association Panel", "GORP", "Biparental"), each = 100))
)

# Create the plot
ggplot(ld_data, aes(x = Genetic_Distance, y = LD, color = Population)) +
  geom_line(size = 1.2) +
  labs(
    title = "LD Decay (R²) vs Genetic Distance",
    x = "Genetic Distance (cM)",
    y = "Linkage Disequilibrium (R²)"
  ) +
  scale_color_manual(values = c("green", "blue", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element
    
    