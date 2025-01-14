theme <- theme(plot.background = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank(),
               axis.line = ggplot2::element_line(linewidth = 0.2),
               plot.title = ggplot2::element_text(hjust = 0.5,
                                                  face = 'bold',
                                                  size = 12),
               axis.text = ggplot2::element_text(size  = 12,
                                                 color = 'black'),
               axis.title = ggplot2::element_text(size  = 12,
                                                  face = 'bold',
                                                  color = 'black'))

# Expects a dataframe with 2 columns:
# traitValA, traitValB (for type "CONTOUR"),
# with a 3rd column (fitness) (for type "SURFACE")
overlayWalkOnLandscape <-function(df,
                                  type="CONTOUR",
                                  fitCalc,
                                  traitAMin=-1,
                                  traitAMax=1,
                                  traitBMin=-1,
                                  traitBMax=1) {
  
  fitness_x = seq(traitAMin,traitAMax, by=0.05)
  fitness_y = seq(traitBMin,traitBMax, by=0.05)
  fitness_z = outer(fitness_x,fitness_y,fitCalc)
  
  if (type == "CONTOUR") {
    fig <- plot_ly() %>%
      layout(xaxis = list(title = "Trait A", constrain = "domain"),
             yaxis = list(title = "Trait B", scaleanchor="x")) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='contour',
        colorscale = "RdBu",
        colorbar=list(title = "w"),
        line = list(color = 'black', width = 1),
        opacity=1) %>%
      add_trace(
        fig,
        df,
        name = n.popSize,
        x = df$traitValA,
        y = df$traitValB,
        type='scatter',
        mode = 'lines',
        line = list(color = 'black', width = 4, dash = 'solid'),
        opacity = 1)
    return (fig)
  } else if (type == "SURFACE"){
    fig <- plot_ly() %>%
      layout(scene = list(xaxis = list(title = "Trait A"),
                          yaxis = list(title = "Trait B"),
                          zaxis = list(title = "w"),
                          aspectmode='cube')) %>%
      add_trace(
        fig,
        df,
        name = n.popSize,
        x = df$traitValA,
        y = df$traitValB,
        z = df$fitness,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        line = list(width = 10)
      ) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='surface',
        colorbar=list(title = "w"),
        colors = "PuBuGn",
        opacity=0.9)
    return (fig)
  } else {
    print("Type Not Supported")
  }
}

# Plot a population on a a 3d fitness chart
# only works for populations with 2 traits
# TODO: figure out how to add colorbar label
plot3dPopulationFitness <- function(pop, fitCalc) {
  popSize <- nInd(pop)
  df <- data.frame(traitA=numeric(popSize),
                   traitB=numeric(popSize),
                   fitness=numeric(popSize))
  for (i in 1:popSize) {
    df$traitA[i] <- pheno(pop)[i,1]
    df$traitB[i] <- pheno(pop)[i,2]
    df$fitness[i] <- fitCalc(df$traitA[i], df$traitB[i])
  }
  fig <- plot_ly()
  
  plot_ly() %>% 
    layout(scene = list(xaxis = list(title = "Trait 1"),
                        yaxis = list(title = "Trait 2"),
                        zaxis = list(title = "w"),
                        aspectmode='cube')) %>%
    add_trace(
      fig,
      df,
      name = popSize,
      x = df$traitA,
      y = df$traitB,
      z = df$fitness,
      type = 'scatter3d',
      mode = 'markers',
      color=df$fitness)
}

# Plots the trait architecture on a genetic basis
# Works for 1 or 2 trait populations
plotTraitArchitecture <- function(pop, methodType="MethodB", fitFunc) {
  geno <- pullQtlGeno(pop,1)
  if (pop@nTraits > 1) {
    for (t in 2:pop@nTraits) {
      geno <- cbind(geno, pullQtlGeno(pop, trait=t))
    }
  }
  cols <- colnames(geno)
  nLoci <- length(cols)
  popSize <- nrow(geno)
  eff_sizes <- data.frame(id = character(nLoci),
                          eff_size=numeric(nLoci))
  for (l in 1:nLoci) {
    id <- cols[l]
    eff_sizes$id[l] <- id
    locus = geno[,l]
    eff_sizes$eff_size[l] <- getEffectSize(locus,
                                           id,
                                           twoTraitFitFunc,
                                           pop,
                                           methodType)
    
  }
  eff_sizes <- eff_sizes[apply(eff_sizes!=0, 1, all),]
  eff_sizes <- na.omit(eff_sizes)
  eff_sizes <- unique(eff_sizes)
  e <- eff_sizes
  g <- ggplot(data=eff_sizes, aes(x=reorder(id, -eff_size), y=eff_size)) +
    geom_bar(stat="identity") +
    labs(x = "Variant Id", y = "Effect Size") +
    theme +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 0.5,
                                     hjust=1,
                                     size = 6,
                                     margin = margin(b = 10)),
          axis.text.y = element_text(margin = margin(l=10, r=10)))
  return (g)
}

# Show the fitness of a population
plot_fitness <- function(df) {
  ggplot(data=df, aes(x=gen, y=fitness)) +
    geom_line() +
    geom_point() +
    labs(x = "Generation", y = "Fitness")
}

# Plot Genetic Values for Two Traits
plot_hist <- function(pop) {
  gv_a = gv(pop)[,1]
  gv_b = gv(pop)[,2]
  idx = c(1:length(gv_a))
  df <- data.frame(gv_a, gv_b)
  df <- melt(as.data.table(df))
  ggplot(df, aes(x=value, color=variable)) + geom_histogram(binwidth=1, position='identity')
}