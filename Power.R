#Load packages
library(ggplot2)
library(dplyr)

#Adapted from Boutchet et.al (2017)
# x-axis: effect size (pqB²)
pqB2 <- c(0.05, 0.10, 0.15, 0.20, 0.25)

# Hypothetical power values for each method & heritability
df <- data.frame(
  pqB2 = rep(pqB2, times = 4),
  Power = c(
    # GWAS h²=0.7
    0.05, 0.15, 0.30, 0.40, 0.45,
    # GORP h²=0.7
    0.10, 0.30, 0.55, 0.70, 0.85,
    # GWAS h²=0.4
    0.02, 0.05, 0.10, 0.15, 0.18,
    # GORP h²=0.4
    0.05, 0.15, 0.25, 0.35, 0.45
  ),
  Method = rep(c("GWAS h2=0.7", "GORP h2=0.7", "GWAS h2=0.4", "GORP h2=0.4"), each = length(pqB2))
)

# Add columns for aesthetic mapping
df <- df %>%
  mutate(
    MethodType = ifelse(grepl("GWAS", Method), "GWAS", "GORP"),
    h2 = ifelse(grepl("0.7", Method), 0.7, 0.4)
  )

# Plot
ggplot(df, aes(x = pqB2, y = Power, color = MethodType, linetype = factor(h2))) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("GWAS" = "black", "GORP" = "green3")) +
  scale_linetype_manual(values = c("0.7" = "solid", "0.4" = "dashed"),
                        labels = c("h²=0.4", "h²=0.7")) +
  labs(x = expression(pqB^2), y = "Power",
       color = "Method", linetype = "Heritability") +
  theme_bw(base_size = 14) +
  ggtitle("N line = 400")

