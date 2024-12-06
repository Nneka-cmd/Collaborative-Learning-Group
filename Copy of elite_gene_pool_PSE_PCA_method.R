

################################################################################
# DO NOT SKIP ANY LINE BELOW!!!!!!!!!!!!!  
# Created functions are not generic, they rely on simulation parameters!!
# Restart R and clean R's memory before starting each simulation scenario
################################################################################

rm(list = ls()) # Clean R's memory by removing all objects

# Set working directory; don't run this, set yours
# setwd("G:/.shortcut-targets-by-id/1-SPOC0ZYfUrZ5CRzW6YgicQR6Eub_HPq/Alex lab folder/AlphaSimR/Safy_data ")
setwd("/Users/awkena/Desktop/Elite_genepools")
# list.files()



{
  # Load packages
  library(AlphaSimR) # for simulation founder pops and breeding schemes
  # library(adegenet) # for dapc
  library(dplyr) # for results summary
  library(ggplot2) # for making pretty plots
  library(ggpubr)
  library(RColorBrewer) # Colors for groups for the plots
  library(gridExtra) # Arranging ggplots into one pdf file
  library(FactoMineR) # PCA for pop structure
  library(factoextra)
  # library(ComplexHeatmap)
  # library(pheatmap)
  library(dendextend)
  library(ggpubr) 
  library(rstatix)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SIMULATION PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Pearl millet 2n = 2x = 14 chromosomes. I suggest the number of segregating 
#' QTLs we simulate should be in multiples of 7 (21 for oligogenic and 735 for
#' polygenic)

#' For  realistic founder pop simulation parameters, the default values for 
#' GENERIC were used. It also makes simulation faster.
#' 
#' The genetic length of each chr = 1 Morgan 
#' Physical length for each chr = 1x10^8 base pairs. 
#' 
#' Sequences for each chromosome were generated using the 
#' Markovian Coalescent Simulator. 
#' 
#' Recombination rate = 1 Morgans / 1x10^8 base pairs = 1x10^-8 per base pair
#'  -- inferred from genome size
#' Mutation rate was set to 2.5x10-8 per base pair. 
#' 
#' Effective population size was set to 100, with linear piecewise increases 
#' to 500 at 100 generations ago, 1,500 at 1,000 generations ago, 6,000
#' at 10,000 generations ago, and 10,000 at 100,000 generations ago.
#' 

{
  # Defining other simulation parameters
  runs <- 10 # number of iterations/ 5 to test code for now
  n.ind <- 200 # number of individuals in founder pop
  nchr <- 7 # number of chromosomes
  n.site <- 1000 # Number of segregation sites per chromosome
  total_site <- nchr * n.site # Total number of loci across chromosomes
  
  varE <- 10 # Error variance for setting simParam
  
  Ne <- 20
  
  # simulation for low (0.2), moderate (0.5) and high (0.8)
  Ha <- 0.8 # Narrow-sense heritability for phenotyping parents and  F5 RILs
  
  trait.mean <- 55 # Mean of FT in founder pop
  
  n.crossesF1 <- 5 # Number of F1 crosses
  
  n.progenyF1 <- 1 # Number of total progenies per cross in the F1
  
  # Total number of progenies at F2 from all F1 crosses
  n.progenyF2 <- 500
  
  # Number of progeny per F1 cross during F2 pop generation
  n.progeny_ssd <- n.progenyF2/(n.progenyF1*n.crossesF1)
  
  # Number of total QTLs to simulate across chromosomes-
  # nQTLs in multiples of 7 except the monogenic scenario
  # oligogenic (21) vs polygenic loci (735)
  # nQTL <- c(rep(1, times = 2), rep(0, times = 5)) # nQTL = 2
  nQTL <- c(rep(1, times = 4), rep(0, times = 3)) # nQTL = 4
  
  # nQTL <- c(rep(3, times = 6), rep(2, times = 1)) # nQTL = 20
  
  # nQTL <- rep(105, times = 7) # nQTL = 735
  
  total_nQTL <- sum(nQTL)
  
  min.t1 <- 52.5 # minimum value for range selection of trait 1
  
  max.t1 <- 57.5 # maximum value for range selection for trait 1
  
  type <- c("cis", "trans") # Run simulation for cis and trans elites
  
  use_seg_site <- TRUE # Use genome-wide marker data for dapc
  
  inbred <- FALSE # simulate founder pop as outbreds
  
  ngroups <- 5 # Positive integer indicating the number of divergent groups
  
  drift_gen <- 4 * Ne
  
  gen <- 50 # Positive integer indicating the number of generations to randomly mate individuals
  
  # Simulation hypothesis type
  hypothesis <- "PS_PCA" # c("SA", "PSE", "PSR")
  
  # Trait simulation using user-inputted additive effects
  value <- (max.t1 - min.t1 + 3)/total_nQTL
  addeff <- rep(value, total_nQTL) # equal add effect values for each QTL locus
  
  
  
  # Load all functions for the simulation into the global environment
  source("elite_genepool_functions.R")
  
}


{
  #~~~~~~~~~~~~~~~ Empty data structures to hold simulation results ~~~~~~~~~~~~~#
  # Selected cis and trans parents info list object
  cis_parents <- trans_parents <- vector(mode = "list", length = runs)
  names(cis_parents) <- names(trans_parents) <- paste0("run", 1:runs)
  
  # Euclidean distance matrix for cis and trans scenarios for each run
  cis_matx <- trans_matx <- vector(mode = "list", length = runs)
  names(cis_matx) <- names(trans_matx) <- paste0("run", 1:runs)
  
  # Pop assignment for parents after pca for each run
  grp_assign <- vector(mode = "list", length = runs)
  names(grp_assign) <- paste0("run", 1:runs)
  
  
  # Create an empty list to hold results for all simulations 
  pheno_ls <- vector(mode = 'list', length = length(type))
  names(pheno_ls) <- c('cis', 'trans')
  
  
  # Create an empty list to hold results for all simulations
  gv_ls <- vector(mode = 'list', length = length(type))
  names(gv_ls) <- c('cis', 'trans')
  
  
  # Empty data frame to hold results for each nQTL simulation 
  pheno_cis <- pheno_trans <- data.frame(matrix(NA, nrow = n.progenyF2, ncol = runs))
  colnames(pheno_cis) <- colnames(pheno_trans) <-  paste0('run', 1:runs)
  
  
  # Empty data frame to hold results for each nQTL simulation for 20 runs
  gv_cis <- gv_trans <- data.frame(matrix(NA, nrow = n.progenyF2, ncol = runs))
  colnames(gv_cis) <- colnames(gv_trans) <-  paste0('run', 1:runs)
  
  # Extract and store additive effects for QTNs for each run 
  # Extract and store cophenetic distances between genome-wide markers and
  # QTN markers
  add_eff <- cor_phen <- vector(mode = "list", length = runs)
  names(add_eff) <- names(add_eff) <- paste0("run", 1:runs)
  
  # Create new directories to save data and results
  dir <- file.path(getwd(), paste0("nRuns_", runs, "|", "herit_", Ha, "|", "Nsite_", total_site, 
                                   "|", "nQTL_", total_nQTL, "|", "F2_pop_size_", 
                                   n.progenyF2, "|", "hypothesis_", hypothesis, "|",
                                   format(Sys.time(), "%F_%H_%M_%S")))
  
  if (!dir.exists(dir)) dir.create(dir)
  
  dir_dapc_plt <- file.path(dir, "DAPC_plots")
  if (!dir.exists(dir_dapc_plt)) dir.create(dir_dapc_plt)
  
  dir_DAPC_data <- file.path(dir, "DAPC_data")
  if (!dir.exists(dir_DAPC_data)) dir.create(dir_DAPC_data)
  
  dir_geno_data <- file.path(dir, "geno_data")
  if (!dir.exists(dir_geno_data)) dir.create(dir_geno_data)
  
  # Define progres bar parameters
  pro_bar <- txtProgressBar(min = 0, max = runs, style = 3, char = '=')
  
  pro_val <- 1:runs
  #names(pro_val) <- nQTL
  pro_val
}

###############################################################################
# Run simulations using defined parameters and loaded functions above
# Do not run this for loop if you have skipped any line of code above!
###############################################################################

for (i in seq_len(runs)) {
  
  # Founder pop simulation will repeat until one with k = 2 - 9 is obtained
  
  repeat {  
    # create founder pop
    founderPop <- runMacs(nInd = n.ind, 
                          nChr = nchr, 
                          segSites = n.site,
                          inbred = inbred, 
                          species = "GENERIC",
                          manualCommand = paste(
                            "1000000000 -t", #Physical lenght 1e8 base pairs
                            2.5/1E8*(4*Ne), #Mutation rate adjusted for Ne
                            "-r",1/1E8*(4*Ne), #Recombination rate adjusted for Ne
                            "-eN",10/(4*Ne), 100/Ne), #Modeling Ne=100 at 10 generations ago
                          manualGenLen = 2)
    
    # Set global simulation parameters and add trait -- additive effects
    SP <- SimParam$new(founderPop)
    
    # Trait simulation using user-inputted additive effects of equal sizes
    mk_names <- get_marker_names(nchr = nchr,
                                 total_nQTL = total_nQTL,
                                 n.site = n.site)
    
    SP <- SP$importTrait(markerNames = mk_names,
                         addEff = addeff,
                         intercept = trait.mean,
                         varE = varE)
    
    # Set global simulation parameters and add trait -- additive effects
    # nQTL <- sample(nQTL, size = nchr) # QTLs per chr
    # SP <- SP$addTraitA(nQtlPerChr = nQTL, mean = trait.mean, var = varE)
    
    # Simulate base pop for selecting parents
    parents <- newPop(founderPop, simParam = SP)
    
    parents <- setPheno(parents, h2 = Ha, simParam = SP)
    
    # Stabilizing selection for acquired trait before pop divergence
    parents <- parents[which(pheno(parents)[,1] >= min.t1 & pheno(parents)[,1] <= max.t1)]
    
    if (parents@nInd > 10) break
    
  }
    
  # Create divergent populations in base pop with optimized selection for 
    # the acquired trait
    parents <- drift_pop(pop = parents,
                         ngroups = ngroups,
                         gen = gen,
                         drift_gen = drift_gen
                         )
    
    # Get sub populations
    div_grp <- as.data.frame(parents$group, stringsAsFactors = TRUE)
    #View(div_grp)
    
    parents <- parents$pop # Pooled populations
    
    
    # Pull genotype data for individuals 
    # Geno data for all seg site
    newpop_geno <- pullSegSiteGeno(parents) |> as.data.frame()
    
    newpop_geno <- cbind(div_grp, newpop_geno)
    
    # The variables id and subpop (index = 1:2) removed
    # before PCA analysis
    pca1 <- PCA(newpop_geno[,-c(1:2)], 
                ncp = 2,
                graph = FALSE)
    
    # Compute hierarchical clustering on principal components
    res.hcpc <- HCPC(pca1, graph = FALSE)
    
    post_pca_grp <- data.frame(id = rownames(res.hcpc$data.clust),
                               grp = res.hcpc$data.clust[,ncol(res.hcpc$data.clust)])
    
  # Define number of pca identified populations
  ngrps <- length(levels(post_pca_grp$grp))
  
  # Get pca individual group assignment and coordinates
  grps <- data.frame(ids = parents@id, pop = post_pca_grp$grp, 
                     gv1 = round(gv(parents)[,1], 2),
                     trait1 = round(pheno(parents)[,1], 2),
                     pca1$ind$coord)
  
  grp_assign[[i]] <- grps
  
  # myCol <- RColorBrewer::brewer.pal(ngrps, "Set1")
  
  # PCA plot using ggplot2
  # Default plot
  plt1 <- fviz_pca_ind(pca1, 
                       pointsize = 2,
                       label="none", 
                       habillage = newpop_geno$subpop) + theme_classic() 
  
  # Make compoplot using ggplot2 
  plt2 <- fviz_dend(res.hcpc, 
                    cex = 0.3, # Label size
                    type = "circular",
                    palette = "jco",
                    rect = TRUE, 
                    rect_fill = TRUE, 
                    rect_border = "jco",           
                    labels_track_height = 0.8)
  
  #' Make scree plot 
  plt3 <- fviz_screeplot(pca1, ncp = 10) + theme_classic() 
  
  # Save all dapc plots in one pdf file
  plts <- list(plt1, plt3, plt2)
  
  ggsave(filename = paste0("plots_", "run_", i, ".pdf"), 
         plot = marrangeGrob(plts, nrow = 1, ncol = 1), 
         device = "pdf", path = dir_dapc_plt, 
         units = "in", width = 10, height = 5)
  
  
  # Save simulated geno data for each run to a folder
  saveRDS(newpop_geno, file = file.path(dir_geno_data,
                                 paste0("run_", i, ".RDS")))
  
  # Nested for loop to perform simulation for cis and trans elites
  for (j in type) {
    
    if (j == "cis") {
      
      # Subset marker data for selected group
      sel_grp <- sample(post_pca_grp$grp, size = 1)
      
      # Use random_sel to identify cis parents
      cis_mat <- rand_sel(type = "cis")
      
      # Get parents for cis elite crosses
      cross_cis <- as.matrix(cis_mat[, 1:2])
      
      cis_parents[[i]] <- cis_mat
      
      # ssd to F5
      pop <- ssd_pop(type = "cis")
      
      # Set phenotype for FT in F5
      pop <- setPheno(pop, h2 = Ha, simParam = SP)
      
      pheno_cis[,i] <- pheno(pop)[,1]
      gv_cis[,i] <- gv(pop)[,1]
      
    } else {
      
      # Use random_sel to identify trans parents
      trans_mat <- rand_sel(type = "trans")
      
      # Get parents for trans elite crosses
      cross_trans <- as.matrix(trans_mat[, 1:2])
      
      trans_parents[[i]] <- trans_mat
      
      # ssd to F5
      pop <- ssd_pop(type = "trans")
      
      # Set phenotype for FT in F5
      pop <- setPheno(pop, h2 = Ha, simParam = SP)
      
      pheno_trans[,i] <- pheno(pop)[,1]
      gv_trans[,i] <- gv(pop)[,1]
      
    }
    
  }
  
  # set progress bar to chromosome indexing
  setTxtProgressBar(pro_bar, value = pro_val[i])
  cat("           Running", runs, "simulations:", "on Run", i)
  
} # End of simulation for loop

################## Tidying and summarizing simulation results #################

# Tidy simulations results for cis and trans as a list 
{
  pheno_ls[['cis']] <- pheno_cis
  pheno_ls[['trans']] <- pheno_trans
  
  gv_ls[['cis']] <- gv_cis
  gv_ls[['trans']] <- gv_trans
  
  #' Convert parent info list to a data frame using the `df_parents()`
  cis_parents <- df_parents(cis_parents) # for cis
  
  trans_parents <- df_parents(trans_parents) # for trans
  
  #View(cis_parents)
  #View(trans_parents)
  
  # Convert pheno_ls list object to a long format data frame object
  pheno_all <- sim_summ(pheno_ls, n.progeny = n.progenyF2, 
                        type = type, runs = runs)
  
  # Convert gv_ls list object to a long format data frame object
  gv_all <- sim_summ(gv_ls, n.progeny = n.progenyF2, 
                     type = type, runs = runs)
  
  #' Combine simulation results for gv and phenotypic values as one big data frame
  #' using the `comb_all()` function
  dat_all <- comb_all(gv = gv_all, pheno = pheno_all)
  
  
  # Get descriptive statistics for each run
  summ <- dat_all |> dplyr::group_by(data_type, type, runs) |> dplyr::summarise('min' = round(min(value), 2),
                                                                                'max' = round(max(value), 2),
                                                                                'mean_run' = round(mean(value), 2),
                                                                                'sd' = round(sd(value), 2),
                                                                                'var' = round(var(value), 2),
                                                                                'selected' = accept(value),
                                                                                .groups = "keep")
  
}

# View(as.data.frame(summ))

{
  myCol <- RColorBrewer::brewer.pal(3, "Set1")
  myCol <- myCol[-3]
  
  # Trait distribution plot using ggplot2 for pheno
  plt_summ_pheno <- ggplot(data = pheno_all, aes(x = value, color = type,
                                                 group = paste(type, runs))) +
    
    geom_density(linewidth = 0.5) + theme_classic() +
    
    scale_color_manual(values = alpha(myCol, 0.5)) + 
    geom_vline(xintercept = c(min.t1, max.t1), linetype = 2) +
    
    # theme_minimal() + facet_wrap(~ runs, ncol = 5) +
    labs(x = 'Trait', y = 'Density',
         title = paste('Progeny trait distribution for phenotypic values:\n', 
                       'heritability =', Ha, "and", "nQTL =", total_nQTL)) +
    theme(axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 16, color = "black"),
          legend.title = element_text(size = 16, face = "bold"))
  
  plt_summ_pheno
  
  
  # Make barplot for variances
  var_plt <- barplot_var()
  var_plt
  
  
  # Combine all into one pdf file
  plt_summs <- list(plt_summ_pheno, var_plt)
  
  ggsave(filename = paste0("progeny_distribution_plots",  ".pdf"), 
         plot = marrangeGrob(plt_summs, nrow = 1, ncol = 1), 
         device = "pdf", path = dir, units = "in",
         width = 6, height = 5)# Save simulated results to a file
}

# Save summarized simulation data for cis and trans 
saveRDS(dat_all, file = file.path(dir, paste0("sim_results", ".RDS")))

# Save parents used for cis and trans crosses
parents_info <- list(cis = cis_parents, trans = trans_parents)
saveRDS(parents_info, file = file.path(dir, paste0("parents_info", ".RDS")))


# Save group assignments after dapc for all runs
# View(grp_assign$run1)
saveRDS(grp_assign, file = file.path(dir, paste0("group_assignment", ".RDS")))

