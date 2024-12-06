#' Identify cis and trans elites for crosses based on the number of
#' shared QTL alleles between genotypes. 
#' @param geno  numeric genotype data with individuals as rows and snps 
#' as column
#' @param type choose either `cis` or `trans` for cis or trans elites, 
#' respectively
#' 
#' @returns a data frame containing the raw IBS between individuals
#' 
#' @details
#' The function calculates the raw IBS among individuals as a metric for shared 
#' QTL alleles between individuals either within a cluster or between clusters.
#' For cis elites, individuals with the highest raw IBS metric within a defined
#' cluster are selected as parents. Conversely, trans elites are those with the
#' lowest raw IBS between all genotypes.
#' #' 
#' 
#' 
ibs_raw <- function(geno, type = c("cis", "trans")){
  
  ## Calculate raw IBS
  geno <- geno/2
  mat <- 1 - cluster::daisy(geno, metric = "gower", warnType = FALSE)
  
  mat <- round(as.matrix(mat), 2)
  
  
  mat[lower.tri(mat, diag = TRUE)] <- NA # subset upper triangle of matrix
  
  # Melt iBS matrix to a data frame
  mat <- cbind(expand.grid(dimnames(mat)), ibs = as.vector(mat))
  
  mat <- na.omit(mat) # Remove NAs
  
  if (type == "cis") {
    # Add the cluster numbers of genotypes
    mat <- merge(mat, grps, by.x = "Var2", 
                 by.y = "ids", sort = FALSE)
    
    mat <- merge(mat, grps, by.x = "Var1", 
                 by.y = "ids", sort = FALSE, suffixes = c(".2", ".1"))
    
    # Re-arrange columns
    mat <- dplyr::relocate(mat, gv1.1, .before = gv1.2)
    mat <- dplyr::relocate(mat, trait1.1, .before = trait1.2)
    
    mat <- mat[order(mat[,3], decreasing = TRUE),]
    
  } else {
    
    # Add the cluster numbers of genotypes
    mat <- merge(mat, grps, by.x = "Var2", 
                 by.y = "ids", sort = FALSE)
    
    mat <- merge(mat, grps, by.x = "Var1", 
                 by.y = "ids", sort = FALSE, suffixes = c(".2", ".1"))
    
    # Re-arrange columns
    mat <- dplyr::relocate(mat, gv1.1, .before = gv1.2)
    mat <- dplyr::relocate(mat, trait1.1, .before = trait1.2)
    
    mat <- mat[order(mat[,3], decreasing = FALSE),]
    
    
  }
  
  return(mat)
  
}


#' Randomly generate marker names for use in importTrait() in AlphaSimR
#' @param nchr A numeric value for the number of chromosomes
#' @param total_nQTL A numeric value for the total number of QTNs
#' @param n.site A numeric value for the number of segregating sites
#' 
#' @returns A character vector of marker names for QTN loci

get_marker_names <- function (nchr, total_nQTL, n.site) { 
  
  if (total_nQTL <= nchr) {
    
    # Character vector of loci ids
    nloci <- as.character(1:n.site)
    
    nqtn_per_chr <- 1 # number of QTN per chromosome
    
    # Sample and sort chromosomes to have QTN
    chr <- sample(as.character(1:nchr), size = if (total_nQTL <= nchr) total_nQTL) |> sort()
    
    # Empty vector to hold locus names of QTN
    mk_names <- vector(mode = "list", length = length(chr))
    names(mk_names) <- chr
    
    # Loop to randomly choose locus names for QTNs
    for (i in seq_len(length(chr))) {
      
      mk_names[[i]] <- paste0(chr[i], "_", sample(nloci, size = nqtn_per_chr))
      
    }
    
  } else if (total_nQTL > nchr && total_nQTL %% nchr == 0) {
    
    nqtn_per_chr <- total_nQTL/nchr # number of QTN per chromosome
    
    # Character vector of loci ids
    nloci <- as.character(1:n.site)
    
    # Sample and sort chromosomes to have QTN
    chr <- as.character(1:nchr) 
    
    # Empty vector to hold locus names of QTN
    mk_names <- vector(mode = "list", length = length(chr))
    names(mk_names) <- chr
    
    # Loop to randomly choose locus names for QTNs
    for (i in seq_len(length(chr))) {
      
      mk_names[[i]] <- sort(paste0(chr[i], "_", sample(nloci, size = nqtn_per_chr)))
      
    }
    
  } else if (total_nQTL > nchr && total_nQTL %% nchr != 0) {
    
    # Minimum number of QTNs per chromosome
    min_nqtn_per_chr <- (total_nQTL - total_nQTL %% nchr)/nchr
    
    # Number of QTNs for last chromosome
    nqtn_last_chr <- min_nqtn_per_chr + total_nQTL %% nchr
    
    nqtn_per_chr <- c(rep(min_nqtn_per_chr, nchr-1), nqtn_last_chr)  # number of QTN per chromosome
    
    # Character vector of loci ids
    nloci <- as.character(1:n.site)
    
    # Sample and sort chromosomes to have QTN
    chr <- as.character(1:nchr) 
    
    # Empty vector to hold locus names of QTN
    mk_names <- vector(mode = "list", length = length(chr))
    names(mk_names) <- chr
    
    # Loop to randomly choose locus names for QTNs
    for (i in seq_len(length(chr))) {
      
      mk_names[[i]] <- sort(paste0(chr[i], "_", sample(nloci, size = nqtn_per_chr[i])))
      
    }
    
  }
  
  return(unname(unlist(mk_names)))
  
}



#' Generate a geometric series additive effects for use in importTrait() in 
#' AlphaSimR randomly or manually
#' 
#' @param method Choose method for generating additive effects: options are 
#' 'random' or 'manual'. The 'random' method uses the rnorm() to sample additive
#' effects from the standard normal distribution. The 'manual' method requires
#' a user-supplied additive effect value.
#' @param total_nQTL A numeric value for the total number of QTNs
#' @param value A numeric value of additive effects to be provided if
#' method = 'manual'.
#' @param scale A logical value. If TRUE, generated additive effects are
#' scaled to Z score.
#' @param scale A logical value: if TRUE, values are scaled to Z scores instead
#' of a geometric series. This is to ensure all additive effect values are 
#' within biological possibilities.
#' @param threshold A numeric value for a threshold additive effect value. Set what is
#' biologically real for the trait.
#' 
#' @returns A numeric vector of geometric sequence of additive effects.


get_add_geom <- function (method = c("random", "manual"), 
                          total_nQTL, 
                          value = NULL,
                          threshold = 10) {
  
  # Get value if method == "random" from the standard normal distribution
  efft <- if (method == "random") {
    
    repeat {
      
      eff <- abs(rnorm(1))
      
      # Value should be within acceptable threshold for trait
      if (max(eff) <= threshold && min(eff) >= (-1*threshold)) {
        
        break
        
      }
    }
    
    eff
    
  } else {
    
    if (!is.null(value)) abs(value) else stop("Provide additive effect value!")
    
  }
  
  # Geometric order of additive effects
  addEffs <- sort(efft^(1/seq(from = 1, to = total_nQTL, by = 1)),
                  decreasing = TRUE)
  
  
  addEffs <- ifelse(addEffs <= median(addEffs), -addEffs, addEffs)
  
  
  return(addEffs)
  
}



get_add_equal <- function(total_nQTL, 
                          value = NULL) {
  
  if (!is.null(value)) abs(value) else stop("Provide additive effect value!")
  
  # Create additive effects of equal sizes
  addEffs <- rep(value, times = total_nQTL)
  
  # Throw in some positive and negative effects
  if (length(total_nQTL) %% 2 == 0) {
    
    
    
  }
  
  return(addEffs)
  
  
}

#' Compute the proportion of acceptable range of FT among progenies
#' @param x is a numeric vector of phenotypic values
#' @param minm is a numeric value of minimum phenotypic value that is acceptable
#' @param maxm is a numeric value of maximun phenotypic value that is acceptable
#' @returns percentage of individuals within the acceptable range of phenotype

accept <- function(x, minm = 52, maxm = 57){
  
  len <- length(x[x >= minm & x <= maxm])
  
  pct <- len/length(x) 
  
  pct <- noquote(scales::percent(pct))
  
  return(pct)
  
}

# accept(pheno$FLO)

#' Summarize simulated results
#' @param x list object from alphasimR simulation for cis and trans elite
#' genotypes
#' @param n.progeny is the number of progenies in x
#' @param type a vector containing cis and trans elements
#' @param runs is the number of replications or runs used in the simulation
#' 
#' @returns a melted data frame object for plotting
#' 
#' 
sim_summ <- function(x, n.progeny, type, runs) {
  
  # Create a vector of nQTLs simulated with each repeated n.progeny times
  type <- rep(type, each = n.progeny)
  
  res <- do.call(rbind, x)
  colnames(res) <- paste0("run", 1:runs)
  
  res <- cbind(type = type, res)
  
  res <- reshape2::melt(res, id.vars = 1,
                        variable.name = "runs", value.name = "value")
  
  res$type <- as.factor(res$type)
  res$runs <- as.factor(res$runs)
  
  return(res)
}


#' Convert list of parent information to a data frame
#' @param x is a list object containing information on parents used in 
#' the crosses
#' 
#' @returns a data frame of parent information
#' 
df_parents <- function(x){
  
  df <- do.call(rbind, x)
  df$Var1 <- as.character(df$Var1)
  df$Var2 <- as.character(df$Var2)
  
  #rownames(df) <- paste0("run_", 1:nrow(df))
  
  return(df)
  
}


# Combine gv and pheno data frames into one data frame
#' @param gv is data frame for gv results for progenies
#' @param pheno is data frame for pheno results for progenies
#' 
#' @returns a data frame containing gv and pheno data.   #### Modify to summarize the variance across runs in another figure
#' 
comb_all <- function(gv, pheno) {
  
  gv <- cbind(data_type = 'gv', gv)
  
  pheno <- cbind(data_type = 'pheno', pheno)
  
  dat_all <- rbind(gv, pheno)
  
  dat_all$data_type <- as.factor(dat_all$data_type)
  
  dat_all$type <- as.factor(dat_all$type)
  
  return(dat_all)
  
}

# Make bar plot for cis and trans variances
# This reduces the number of times to two instead of running multiple lines

barplot_var <- function () {
  
  sum_pheno <- as.data.frame(dplyr::filter(summ, data_type == 'pheno'))
  names(sum_pheno)[8] <- "var.a"
  ### Calculate the total of variance for all the mean
  ### mean variance plot 
  #Extract cis elite crosses. 
  extract_cis <- dplyr::filter(sum_pheno, type == 'cis')
  #Extract trans elite crosses. 
  extract_trans <- dplyr::filter(sum_pheno, type == 'trans')
  
  tt <- t_test(data = sum_pheno, var.a ~ type, alternative = "less",
               p.adjust.method = 'bonferroni',
               ref.group = "cis") |> add_significance() |> add_xy_position(x = "type")
  
  #Mean 
  mean_var <- c(cis = mean(extract_cis$var), trans = mean(extract_trans$var))
  
  #standard deviation and error
  sd_var <- c(cis = sd(extract_cis$var) , trans = sd(extract_trans$var))
  
  se_var <- sd_var/sqrt(runs)
  
  # Put all summary data for variance together into a data frame
  barplot_dat <- data.frame(Elite = c("cis", "trans"),  Variance = mean_var,
                            sd = sd_var, se = se_var)
  
  # Barplot 
  P <- ggplot2::ggplot(barplot_dat, ggplot2::aes(x = Elite, y = Variance, color = Elite)) + 
    
    ggplot2::geom_bar(stat = "identity", fill = "white", linewidth = 1) +
    
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Variance, ymax = Variance + sd), 
                           width = 0.2,  position = ggplot2::position_dodge(0.05)) +
    
    ggplot2::scale_color_manual(values = myCol) + ylim(c(0, max(barplot_dat$Variance) + max(barplot_dat$sd) + 0.08)) +
    
    ggplot2::labs(x = 'Elite type', y = 'Variance', 
                  title = paste('Variance in cis vs. trans-elites: \n','runs =', runs)) +
    
    ggplot2::theme_classic() +  theme(axis.text = element_text(size = 16, color = "black"),
                                      axis.title = element_text(size = 16, face = "bold"),
                                      legend.text = element_text(size = 16, color = "black"),
                                      legend.title = element_text(size = 16, face = "bold"))
  
  P + stat_pvalue_manual(tt, 
                         label = "p.signif", 
                         size = 10, 
                         y.position = max(barplot_dat$Variance) + max(barplot_dat$sd) + 0.05)
  
  
}

#' Compute the number of loci with different genotypes for cis and 
#' trans parents
#' @param type Set to cis or trans accordingly.
#' @param nselect An integer indicating the number of rows of matrix to select.
#' @returns A numeric vector of the number of loci with different genotypes.

allele_diff <- function(type = c("cis", "trans"),
                        nselect) {
  
  if (type == "cis") {
    
    dat <- cis_mat
    
  } else (dat <- trans_mat)
  
  par1 <- dat$Var1[1:nselect]
  par2 <- dat$Var2[1:nselect]
  
  diff <- vector(mode = "numeric", length = nselect)
  
  for (n in seq_len(nselect)) {
    
    par1_geno <- as.vector(qtl_geno[par1[n],])
    par2_geno <- as.vector(qtl_geno[par2[n],])
    
    comp <- par1_geno == par2_geno
    
    diff[n] <- length(comp[comp == FALSE])
    
  }
  
  return(diff)
}



# Perform single seed descent to F5
#' @returns An F5 AlphasimR pop object
ssd_pop <- function(type = c("cis", "trans")) {
  
  if(type == "cis") {
# Make crosses using selected pair of cis elite parents
F1 <- makeCross(pop = parents, crossPlan = cross_cis, 
                nProgeny = n.progenyF1, simParam = SP)
  } else {
  
    # Make crosses using selected pair of cis elite parents
    F1 <- makeCross(pop = parents, crossPlan = cross_trans, 
                    nProgeny = n.progenyF1, simParam = SP)
}
# Self F1 to F2 
pop <- self(F1, nProgeny = n.progeny_ssd, simParam = SP)


# Self F2 to F3 by SSD
pop <- self(pop, nProgeny = 1, simParam = SP)

# Self F3 to F4 by SSD
pop <- self(pop, nProgeny = 1, simParam = SP)

# Self F4 to F5 by SSD
pop <- self(pop, nProgeny = 1, simParam = SP)

return(pop)

}


#' Perform recurrent selection on cis and trans F5 progenies on desired traits
#' @param SI proportion of individuals to select
#' @param n.gen Number of cycles/generations of recurrent selection to run
#' @param n.progenyF2 Total number of progenies at F2 or F5


recurrent_sel <- function(SI, n.gen, pop) {
  
  # nSelected <- floor(length(pheno(pop)[,1]) * SI)# # number of individuals selected in the F5 progenies based on the selection intensity 
  
  # Function to calculate mean and variance on pop object (x)
  pop_sum <- function(x) {
    
    dat <- c(meanP(x), meanG(x), diag(varG(x)), diag(varP(x)))
    dat
  }
  
  
  
  #Store data for each generation 
  newPop_data <- data.frame(matrix(NA, nrow = n.gen+1, ncol = 9))
  
  colnames(newPop_data) <- c("cycle", "mean.trait1", "mean.trait2", 
                             "mean.gv1", "mean.gv2", 
                             "VarG1", "VarG2", "VarP1", "VarP2")
  
  rownames(newPop_data) <- paste0("cycle_", 0:n.gen)
  
  newPop_data[,1] <- 0:n.gen
  
  newPop_data[1,2:ncol(newPop_data)] <- pop_sum(parents)
  
  # Empty list objects to hold results for each run of simulation for phenotypic values 
  # pheno_1 <- pheno_2 <- vector(mode = "list", length = n.gen)
  # names(pheno_1) <-  names(pheno_2) <- paste0('cycle', 1:n.gen)
  
  # Empty data frame to hold results for each run of simulation for phenotypic values 
  pheno_1 <- pheno_2 <- data.frame(matrix(NA, nrow = n.progenyF2, ncol = n.gen))
  colnames(pheno_1) <-  colnames(pheno_2) <- paste0('cycle', 1:n.gen)
  
  
  for (g in seq_len(n.gen)) {
    
    # Range selection for trait 1 before recurrent selection
    pop <- pop[which(pheno(pop)[,1] >= min.t1 & pheno(pop)[,1] <= max.t1)]
    
    nSelected <- floor(length(pheno(pop)[,1]) * SI)# # number of individuals selected in the F5 progenies based on the selection intensity
    
    newParents <- selectInd(pop = pop,
                            nInd = nSelected,
                            trait = 2,
                            use = "pheno",
                            selectTop = TRUE,
                            returnPop = TRUE)
    
    if (newParents@nInd < 2) {
      
      break
      
    }
    
    # Select parents that fall within the acceptable threshold for acquired traits
    # before random cross
    # pop1 <- sel_range(newParents = newParents, min = 50, max = 55)
    
    # Random mating of new population  after selection
    new_F1 <- randCross(pop = newParents, nCrosses = n.crossesF1,
                        nProgeny = n.progenyF1, 
                        simParam = SP)
    
    # Self F1 to F2 
    pop <- self(new_F1, nProgeny = n.progeny_ssd, simParam = SP)
    
    # Self F2 to F3 by SSD
    pop <- self(pop, nProgeny = 1, simParam = SP)
    
    # Self F3 to F4 by SSD
    pop <- self(pop, nProgeny = 1, simParam = SP)
    
    # Self F4 to F5 by SSD
    pop <- self(pop, nProgeny = 1, simParam = SP)
    
    # Set phenotype for FT in F5
    pop <- setPheno(pop, h2 = c(H1, H2), simParam = SP) ### Do we want to have the same heritability for both traits????
    
    # Get mean and variance for each cycle
    newPop_data[(g+1), 2:ncol(newPop_data)] <- pop_sum(pop)
    
    # pheno_1[[g]] <- pheno(pop)[,1]
    # pheno_2[[g]] <- pheno(pop)[,2]
    
    pheno_1[,g] <- pheno(pop)[,1]
    pheno_2[,g] <- pheno(pop)[,2]
    
  }
  
  res <- list(mean_var = newPop_data, trait1 = pheno_1, trait2 = pheno_2)
  
  return(res)
  
}


# Perform dapc using denovo pop clusters
#' @param geno a data matrix of pulled genotypes at all segregating sites 
#' simulated for founder population
#' 
#' @returns all dapc output using the `adegenet` package
#' @details
#' The population divergence grouping is used as the prior groups before dapc.
#' 
#' 
dapc_denovo <- function(geno) {
  
  # # find denovo groups using kmeans before dapc
  find.grp <- adegenet::find.clusters(x = geno, max.n.clust = 10, n.pca = 500,
                                     choose.n.clust = FALSE, criterion = "diffNgroup")

  # Run first dapc using divergence groups in founder pop
  dapc.groups <- adegenet::dapc(geno, var.contrib = TRUE, scale = FALSE, 
                                n.pca = n_pca, grp = find.grp$grp, n.da = 5)
  
  # Find the optimal number of pcs to use
  optim <- adegenet::optim.a.score(x = dapc.groups, 
                                   n.pca = 1:ncol(dapc.groups$tab),
                                   plot = FALSE, n.sim = 100)
  
  # Re-run second dapc using optimal number of pcs
  optim.dapc <- adegenet::dapc(geno, var.contrib = TRUE, scale = FALSE, 
                               n.pca = optim$best, grp = find.grp$grp, n.da = 5)
  
  return(optim.dapc)
}


# The slot assign.per.pop indicates the proportions of successful 
# reassignment (based on the discriminant functions) of individuals 
# to their original clusters. Large values indicate
# clear-cut clusters, while low values suggest admixed groups

# Plot accumulated variance of pcs used in dapc
scree_plot <- function(){
  
  len <- length(optim.dapc$pca.eig)
  temp <- data.frame(ids = 1:len, 
                     eig.val = 100 * cumsum(optim.dapc$pca.eig)/sum(optim.dapc$pca.eig),
                     group = rep(c("black", "lightgrey"),
                                 c(optim.dapc$n.pca, (len - optim.dapc$n.pca))))
  
  dd <- ggplot2::ggplot(temp, ggplot2::aes(x = ids, y = eig.val, fill = group)) +
    
    ggplot2::geom_bar(stat = "identity", linewidth = 1, position = "stack") +
    ggplot2::geom_text(ggplot2::aes( x = 25, y = 90,
                                     label = paste("PC variance used =", round(optim.dapc$var, 2))))+
    
    ggplot2::labs(x = "PCA axis", y = "Acumulated variance (%)") +
    ggplot2::xlim(c(0, len+1)) +  ggplot2::ylim(c(0, 100)) + ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = c("black", "lightgrey")) +
    ggplot2::theme(legend.position = "none")
  
  dd
  
}

#' Identify cis and trans elites for crosses based on the Euclidean
#' distance between genotypes within a cluster (for cis) or genotypes between
#' clusters (for trans). 
#' @param geno  numeric genotype data with individuals as rows and snps 
#' as column
#' @param type choose either `cis` or `trans` for cis or trans elites, 
#' respectively
#' 
#' @returns a data frame containing the Euclidean distance between individuals
#' 
#' @details
#' The function calculates the Euclidean distance among individuals as a metric for
#' genetic relatedness within a cluster or between clusters.
#' For cis elites, individuals with the lowest distance are selected as parents. 
#' Conversely, trans elites are those with the highest distance between all genotypes.
#' 
#' 
#' 
euclid_dist <- function (geno, 
                        type = c("cis", "trans"), 
                        method = c("euclid", "cor")) {
  
  if (type == "cis") {
    
    geno <- geno[grps$pop == pop_cis,]
    
  } else {
    
    geno <- geno
    
  }
  
  if (method == "euclid") {
    
    ## Calculate Euclidean distance
    mat <- round(as.matrix(dist(geno)), 2)
    
  } else {
    geno <- t(geno)
    mat <- as.dist((1 - cor(geno))/2)
    mat <- round(as.matrix(mat), 2)
  }
  
  mat[lower.tri(mat, diag = TRUE)] <- NA # subset upper triangle of matrix
  
  # Melt Euclidean distance matrix to a data frame
  mat <- cbind(expand.grid(dimnames(mat)), dis = as.vector(mat))
  
  mat <- na.omit(mat) # Remove NAs
  
  # Add the cluster numbers of genotypes
  mat <- merge(mat, grps[,c(1:4)], by.x = "Var2",
               by.y = "ids", sort = FALSE)
  
  mat <- merge(mat, grps[,c(1:4)], by.x = "Var1",
               by.y = "ids", sort = FALSE, suffixes = c(".2", ".1"))
  
  # Re-arrange columns
  mat <- dplyr::relocate(mat, pop.1, .before = pop.2)
  mat <- dplyr::relocate(mat, gv1.1, .before = gv1.2)
  mat <- dplyr::relocate(mat, trait1.1, .before = trait1.2)
  
  if (type == "cis") {
    
    mat <- mat[order(mat[,3], decreasing = FALSE),]
    
  } else {
    
    mat <- mat[order(mat[,3], decreasing = TRUE),]
    
  }
  
  return(mat)
  
}


#' Identify cis and trans parents randomly for crossing
#' @param type choose either `cis` or `trans` to select parents for cis or trans
#' elites, respectively.
#' 
#' @returns a data frame containing the selected cis or trans elite parents
#' 
#' @details
#' For cis elite parents, individuals are selected randomly from a given cluster.
#' For trans elite parents, individuals are selected randomly from two different
#' clusters with the highest inter-cluster distance.
#' 
#' 

rand_sel <- function(type = c("cis", "trans")) {
  
  # Create matrix to hold selected parents
  mat <- as.data.frame(matrix(NA, nrow = n.crossesF1, ncol = 2))
  colnames(mat) <- c("Var1", "Var2")
  
  if (type == "cis") {
    
    # Select groups for cis parents
    cis_grp <- grps$ids[grps$pop == sel_grp]
    
    for (m in seq_len(n.crossesF1)) {
      
      mat[m,] <- sample(cis_grp, size = 2, replace = FALSE)
      
      
    }
    
    mat <- cbind(mat, pop = sel_grp) # add pop assignment
    
    # Add the cluster numbers of genotypes
    mat <- merge(mat, grps[,c(1, 3:4)], by.x = "Var2",
                 by.y = "ids", sort = FALSE)
    
    mat <- merge(mat, grps[,c(1, 3:4)], by.x = "Var1",
                 by.y = "ids", sort = FALSE, suffixes = c(".2", ".1"))
    
    # Re-arrange columns
    mat <- dplyr::relocate(mat, gv1.1, .before = gv1.2)
    mat <- dplyr::relocate(mat, trait1.1, .before = trait1.2)
    
    
  } else {
    
    # Inter-cluster distance using grouping coordinates from dapc -- faster
    dis <- round(as.matrix(dist(pca1$ind$coord)), 4)
    
    dis[lower.tri(dis, diag = TRUE)] <- NA # subset upper triangle of matrix
    
    # Melt inter-cluster distance matrix to a data frame
    dis <- cbind(expand.grid(dimnames(dis)), pop_dis = as.vector(dis))
    
    dis <- na.omit(dis) # Remove NAs
    
    sel_pops <- dis[which.max(dis$pop_dis),] 
    
    trans_grp1 <- grps$pop[grps$ids == sel_pops[1,1]]
    trans_grp2 <- grps$pop[grps$ids == sel_pops[1,2]]
    
    
    for (m in seq_len(n.crossesF1)) {
      
      mat[m, 1] <- sample(grps$ids[grps$pop == trans_grp1], size = 1)
      mat[m, 2] <- sample(grps$ids[grps$pop == trans_grp2], size = 1)
      
    }
    
    # Add the cluster numbers of genotypes
    mat <- merge(mat, grps[,c(1:4)], by.x = "Var2",
                 by.y = "ids", sort = FALSE)
    
    mat <- merge(mat, grps[,c(1:4)], by.x = "Var1",
                 by.y = "ids", sort = FALSE, suffixes = c(".2", ".1"))
    
    # Re-arrange columns
    mat <- dplyr::relocate(mat, pop.1, .before = pop.2)
    mat <- dplyr::relocate(mat, gv1.1, .before = gv1.2)
    mat <- dplyr::relocate(mat, trait1.1, .before = trait1.2)
    
  }
  
  return(mat)
  
}

#' Create divergent populations from a founder pop.
#' @param pop AlphasimR pop object from founder pop.
#' @param gen  Positive integer indicating the number of generations to randomly
#' mate individuals.
#' @param ngroups Positive integer indicating the number of divergent groups.
#' @param N An integer indicating the number of individuals to sample to create
#' a founder effect
#' @details
#' The function creates divergent sub-populations based on a defined number of 
#' sub-groups and number of generations of random mating with optimized selection.
#' Optimized selection is implemented by a predefined range of values within the 
#' acceptable range of values for the acquired trait.
#' 

diverge_pop <- function(pop, 
                        gen = 50, 
                        ngroups = 5,
                        N = 10){
  
  # Parent pop trim parameter if pop size is not a multiple of ngroups
  trim_par <- pop@nInd %% ngroups  
  
  if (trim_par > 0) {
    
    st_vec <- pop@nInd - trim_par + 1 # start integer
    pop <- pop[-c(st_vec:pop@nInd)]
    
    
  } else if (trim_par == 0) {
    
    pop <- pop
    
  }
  
  chunk_size <- pop@nInd/ngroups # Positive integer indicating the number of individuals in each sub-group
  
  # Empty list object to hold results for each sub-group
  rand_mate <- div_grp <- vector(mode = "list", length = ngroups)
  
  # Get IDs of individuals in pop object and reshuffle
  iids <- sample(pop@id, size = pop@nInd)
  
  # Create sub-groups from pop object by randomly selecting individuals 
  # without replacement
  splits <- split(iids, rep(1:ngroups, each = chunk_size))
  
  for (i in seq_len(ngroups)) {
    
    # Get each sub-group by sampling individuals without replacement
    # This is to create a founder effect
    grp <- pop[splits[[i]]] 
    
    nn <- 0  # iteration counter
    
    
    repeat {
      
      nn <- nn + 1 # iteration counter
      
      # Random mating in each sub-group with selection
      # Random mating in sub-group
      grp <- randCross(grp, nCrosses = 100, simParam = SP, nProgeny = 1)
      
      grp <- setPheno(grp, h2 = Ha, simParam = SP)
      
      # Range selection for trait 1 before cis/trans parent selection
      grp <- grp[which(pheno(grp)[,1] >= min.t1 & pheno(grp)[,1] <= max.t1)]
      
      # break statement to exit repeat loop
      if (nn == gen) {
        
        break
        
      }
      
    }
    
    # Fill rand_mate object with results
    
    grp <- grp[sample(grp@id, size = chunk_size)]
    
    div_grp[[i]] <- cbind(id = grp@id, subpop = i)
    
    # Randomly select individuals based on chunk size
    rand_mate[[i]] <- grp
  }
  
  # Merge all sub-groups again after random mating and selfing
  newPop <- mergePops(rand_mate)
  
  div_grp_comb <- do.call(rbind, div_grp)
  
  res <- list(pop = newPop, group = div_grp_comb)
  
  return(res)
  
}


#' Create divergent populations from a founder pop.
#' @param pop AlphasimR pop object from founder pop.
#' @param gen  Positive integer indicating the number of generations to randomly
#' mate individuals in the pop object.
#' @param ngroups Positive integer indicating the number of divergent groups.
#' @param drift_gen An integer indicating the number of generations to drift 
#' subpopulations before merging them.
#' @details
#' The function creates divergent sub-populations based on a defined number of 
#' sub-groups and number of generations of random mating with optimized selection.
#' Optimized selection is implemented by a predefined range of values within the 
#' acceptable range of values for the acquired trait.
#' 

drift_pop <- function(pop, 
                      gen = 50,
                      ngroups,
                      drift_gen){
  
  # Empty list object to hold results for each sub-group
  rand_mate <- div_grp <- vector(mode = "list", length = ngroups)
  
  nn <- 0  # iteration counter
  
  # Random mating within the founder before split with stabilizing selection
  repeat {
    
    nn <- nn + 1 # iteration counter
    
    # Random mating in each sub-group with selection
    # Random mating in sub-group
    grp <- randCross(pop, nCrosses = n.ind, simParam = SP, nProgeny = 1)
    
    grp <- setPheno(grp, h2 = Ha, simParam = SP)
    
    # Stabilizing selection for trait 1 before cis/trans parent selection
    grp <- grp[which(pheno(grp)[,1] >= min.t1 & pheno(grp)[,1] <= max.t1)]
    
    # break statement to exit repeat loop
    if (nn == gen) {
      
      break
      
    }
    
  }
  
  # Split and drifting of the subpopulations
  # drift_gen <- 4 * Ne
  for (i in seq_len(ngroups)) {
    
    subpop <- grp[sample(grp@id, size = Ne, replace = TRUE)]
    
    
    dd <- 0 # iteration counter
    
    repeat {
      
      dd <- dd + 1
      subpop <- subpop[sample(subpop@id, size = Ne, replace = TRUE)]
      subpop <- randCross(subpop, nCrosses = 50, simParam = SP, nProgeny = 1)
      
      subpop <- setPheno(subpop, h2 = Ha, simParam = SP)
      
      # Stabilizing selection for trait 1 before cis/trans parent selection
      subpop <- subpop[which(pheno(subpop)[,1] >= min.t1 & pheno(subpop)[,1] <= max.t1)]
      
      # break statement to exit repeat loop
      if (dd == drift_gen) {
        
        break
        
      }
      
    }
    
    rand_mate[[i]] <- subpop
    
    # Get subpop number 
    div_grp[[i]] <- cbind(id = subpop@id, subpop = i)
  }
  
  
  # Merge all sub-groups again after random mating and selfing
  newPop <- mergePops(rand_mate)
  
  div_grp_comb <- do.call(rbind, div_grp)
  
  res <- list(pop = newPop, group = div_grp_comb)
  
  return(res)
  
}



#' Get marker positions and chromosome numbers.
#' @param x Pulled genome-wide or QTL marker loci from AlphaSimR pop object.
#' @returns A data frame of marker loci and their chromosome positions
get_mks <- function(x) {
  
  # Get marker names from imported data
  markers <- colnames(x)
  
  # Parse marker names to extract chromosome numbers and physical positions
  chr <- t(as.data.frame(strsplit(markers, "_"))) 
  chr <- as.data.frame(chr, row.names = markers)
  colnames(chr) <- c('chr', 'pos')
  
  chr$pos <- as.numeric(chr$pos) # Convert marker positions to numeric values
  
  return(chr)
  
}

# Does clustering based on genome-wide marker data depict cis/trans
co_phenetic <- function() {
  
  cor_phe <- numeric()
  
  for (z in grps$pop) {
    
    qtl_grp <- qtl_geno[grps$pop == z,]
    geno_grp <- geno[grps$pop == z,]
    
    if (method == "cor") {
      gg1 <- t(qtl_grp) 
      gg2 <- t(geno_grp)
      qtl_dis <- as.dist((1 - cor(gg1))/2)
      geno_dis <- as.dist((1 - cor(gg2))/2)
      
      hc1_qtl <- hclust(qtl_dis, method = "ward.D2")
      hc1_geno <- hclust(geno_dis, method = "ward.D2")
      
    } else {
      
      hc1_qtl <- hclust(dist(qtl_grp), method = "ward.D2")
      hc1_geno <- hclust(dist(geno_grp), method = "ward.D2")
      
    }
    
    # Co-phenetic correlation coefficient
    D1 <- as.vector(as.dist(cophenetic(hc1_qtl)))
    D2 <- as.vector(as.dist(cophenetic(hc1_geno)))
    
    cor_phe[z] <- cor(D1, D2)
    
  }
  
  return(cor_phe)
}
