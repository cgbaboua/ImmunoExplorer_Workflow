# Functions for WGCNA on CSU Dataset (Author : Cassandra GBABOUA)

# Function to Pick a threshold and plot scale independence and mean connectivity graphs
pickThresholdAndPlot <- function(netwk_type,verbose,return_data=FALSE) {
  # Call the network topology analysis function
  sft = pickSoftThreshold(input_mat, powerVector = powers, networkType = netwk_type, verbose = verbose)
  sft_data <- sft$fitIndices
  # Afficher les deux graphiques côte à côte
  #grid.arrange(plotScaleIndependence(sft_data),plotMeanConnectivity(sft_data),nrow=2)
  return(sft_data)
}

# Function to create the scale independence graph
plotScaleIndependence <- function(sft_data) {
  ggplot(sft_data, aes(x = powers, y = SFT.R.sq)) + geom_point() + geom_line() +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
    labs(x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit,\nsigned R^2", title = "Scale Independence") +
    theme_minimal()
}

# Function to create the mean connectivity graph
plotMeanConnectivity <- function(sft_data) {
  ggplot(sft_data, aes(x = powers, y = mean.k.)) + geom_point() + geom_line() +
    labs(x = "Soft Threshold (power)", y = "Mean Connectivity", title = "Mean Connectivity") +
    theme_minimal()
}

# Function to arrange and plot connectivity and scale independance graph
# Afficher les deux graphiques côte à côte
PlotsAndArrange <- function(sft_data){
  grid.arrange(plotScaleIndependence(sft_data),plotMeanConnectivity(sft_data),nrow=2)
}

# Function to determine the optimal power for WGCNA
determineOptimalPower <- function(df, r2_threshold = 0.8) {
  # Ensure data is sorted by power
  df <- df[df$Power >= 1 & df$Power <= 30, ]
  df <- df[order(df$Power),]
  # Find the smallest power where SFT.R.sq exceeds the threshold
  indices_above_threshold <- which(df$SFT.R.sq >= r2_threshold)
  
  if (length(indices_above_threshold) > 0) {
    # If the threshold is met, use the first power that exceeds the threshold
    optimal_power <- df$Power[indices_above_threshold[1]]
  } else {
    # If the threshold is not met, use the power with the highest R²
    index_of_max_r2 <- which.max(df$SFT.R.sq)
    optimal_power <- df$Power[index_of_max_r2]
  }
  
  # Return the optimal power
  return(optimal_power)
}

# Function to construct network and detect modules
constructNetworkAndDetectModules <- function(netwk_type,picked_power,verbose) {
  picked_power = picked_power
  temp_cor <- cor       
  cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)
  netwk <- blockwiseModules(input_mat, # <= input here
                            
                            # == Adjacency Function ==
                            power = picked_power,                # <= power here
                            networkType = netwk_type,
                            
                            # == Tree and Block Options ==
                            deepSplit = 2, # Détermine la profondeur de la division hiérarchique lors de l'identification des modules. Une valeur plus élevée conduit à une division plus fine des modules.
                            pamRespectsDendro = F, # Indique si l'assignation des modules doit respecter strictement la structure dendrologique. "F" (False) permet plus de flexibilité dans l'assignation des gènes aux modules.
                            # detectCutHeight = 0.75,
                            minModuleSize = 30, # Fixe la taille minimale pour qu'un groupe de gènes soit considéré comme un module. Cela évite de considérer de très petits groupes de gènes comme des modules distincts.
                            
                            maxBlockSize = 4000, # Définit la taille maximale des blocs de données à analyser. Cela est crucial pour gérer efficacement la mémoire et la performance lors de l'analyse de grands ensembles de données.
                            
                            # == Module Adjustments ==
                            reassignThreshold = 0, # Détermine le seuil pour réaffecter les gènes à différents modules lors de l'optimisation du réseau. Une valeur basse implique une réaffectation plus fréquente.
                            mergeCutHeight = 0.25, # Définit le seuil pour la fusion des modules similaires, basé sur la hauteur de coupe dans l'arbre de clustering. Une valeur plus basse conduit à une fusion plus agressive des modules.
                            
                            # == TOM == Archive the run results in TOM file (saves time)
                            saveTOMs = T, # Indique si la Topological Overlap Matrix doit être sauvegardée, ce qui peut être utile pour des analyses ultérieures ou pour économiser du temps de calcul.
                            saveTOMFileBase = "ER", # Spécifie le préfixe de base pour les fichiers de sauvegarde TOM.
                            
                            # == Output Options
                            numericLabels = T, # Détermine si les modules sont étiquetés numériquement ou non. "T" (True) attribue des étiquettes numériques aux modules.
                            verbose = verbose)
  
  cor <- temp_cor     # Return cor function to original namespace
  return(netwk)
}

# Function to merged colors for module visualization
calculateMergedColors <- function(netwk) {
  return(labels2colors(netwk$colors))
}

# Function to visualize modules
visualizeModules <- function(netwk) {
  # Convert labels to colors for plotting
  mergedColors = calculateMergedColors(netwk)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    netwk$dendrograms[[1]], 
    mergedColors[netwk$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05 )
  # netwk$colors[netwk$blockGenes[[1]]]
  # table(netwk$colors)
}

# Function to relates modules to traits and returns module dataframe
relateModulesToTraits <- function(netwk) {
  module_df <- data.frame(
    gene_id = names(netwk$colors),
    colors = calculateMergedColors(netwk)
  )
  
  #write_delim(module_df,
  #file = "gene_modules.txt",
  #delim = "\t")
  
  return(module_df)
}

# Function to check the list of genes by modules
listOfModules <- function(module_df){
  list_of_modules <- split(module_df$gene_id, module_df$colors)
  return(list_of_modules)
}

# Function to merge FC (Fold Change) data with module data
mergeFCWithModuleData <- function(module_df, DEA_FC) {
  module_df_with_fc <- merge(module_df, DEA_FC, by.x = "gene_id", by.y = "Genes")
  module_df_with_fc <- module_df_with_fc[order(module_df_with_fc$colors), ]
  rownames(module_df_with_fc) <- NULL
  return(module_df_with_fc)
}

# Function to summarize module data
summarizeModuleData <- function(module_df) {
  moduleSummary <- setNames(data.frame(table(module_df[, "colors"])), c("ModuleColor", "NbOfGenes"))
  return(moduleSummary)
}

# Function to Plot module-trait relationships
plotModuleTraitRelationships <- function(input_mat, netwk, sample_status,x,y) {
  mergedColors = calculateMergedColors(netwk)
  # Get Module Eigengenes per cluster
  MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
  # Reorder modules so similar modules are next to each other
  if (ncol(MEs0) > 1) {
    MEs0 <- orderMEs(MEs0)
    module_order <- names(MEs0) %>% gsub("ME", "", .)
  } else {
    module_order <- names(MEs0) %>% gsub("ME", "", .)
  }
  
  # Add treatment names
  MEs0$treatment = row.names(MEs0)
  # Add disease_state informations
  MEs0 <- merge(MEs0, sample_status, by.x = x, by.y = y, all.x = TRUE)
  # tidy & plot data
    mME <- MEs0 %>%
      pivot_longer(-c(treatment, DiseaseState)) %>%
      mutate(
        name = gsub("ME", "", name),
        name = factor(name, levels = module_order)
      )
  
  plot <- mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
    geom_tile() +
    facet_grid(~DiseaseState, scales = "free_x", space = "free_x") +
    theme_bw()  +
    scale_fill_gradient2(
      low = "blue", 
      high = "red", 
      mid = "white", 
      midpoint = 0, 
      limit = c(-1,1)) +
    theme(axis.text.x = element_text(angle=90)) +
    labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
  return(list(MEs0 =MEs0,plot= plot))
}

# Function to find the modules where each gene appears
group_genes_by_module <- function(genes_of_interest, list_of_modules) {
  modules <- names(list_of_modules)
  # Create an empty data frame to store the results
  result_df <- data.frame(module = character(), genes = character(), stringsAsFactors = FALSE)
  for (module in modules) {
    # Find the corresponding genes in the current module
    matching_genes <- genes_of_interest[genes_of_interest %in% list_of_modules[[module]]]
    # If genes are found, add them to the data frame
    if (length(matching_genes) > 0) {
      result_df <- rbind(result_df, data.frame(module = module, genes = paste(matching_genes, collapse = ",")))
    }
  }
  # Add a row for genes not found, if necessary
  not_found_genes <- genes_of_interest[!genes_of_interest %in% unlist(list_of_modules)]
  if (length(not_found_genes) > 0) {
    result_df <- rbind(result_df, data.frame(module = NA, genes = paste(not_found_genes, collapse = ",")))
  }
  return(result_df)
}

# Function to examine expression profiles
examineExpressionProfiles <- function(submod,module_df){
  
  subexpr = exprs_data[submod$gene_id,]
  
  submod_df = data.frame(subexpr) %>%
    mutate(
      gene_id = row.names(.)
    ) %>%
    pivot_longer(-gene_id) %>%
    mutate(
      module = module_df[gene_id,]$colors
    )
  
  submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
    geom_line(aes(color = module),
              alpha = 0.2) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    facet_grid(rows = vars(module)) +
    labs(x = "treatment",
         y = "normalized expression")
}

# Function to calculate TOM and display edge list
calculateTOMandEdgeList <- function(module_df, picked_power, modules_of_interest){
  # Filter genes of interest based on selected modules
  genes_of_interest = module_df %>%
    dplyr::filter(colors %in% modules_of_interest)
  
  # Extract expression data for genes of interest
  expr_of_interest = exprs_data[genes_of_interest$gene_id,]
  
  # Only recalculate TOM for modules of interest (faster, though there's some discussion online if this might be slightly off)
  TOM = TOMsimilarityFromExpr(t(expr_of_interest), power = picked_power)
  
  # Add gene names to rows and columns
  row.names(TOM) = row.names(expr_of_interest)
  colnames(TOM) = row.names(expr_of_interest)
  
  # Generate the edge list from TOM
  edge_list = data.frame(TOM) %>%
    dplyr::mutate(gene1 = row.names(.)) %>%
    tidyr::pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
    dplyr::filter(gene1 != gene2) %>%
    dplyr::mutate(module1 = module_df[gene1, "colors"], 
                  module2 = module_df[gene2, "colors"]) %>%
    dplyr::distinct()
  
  # Display edge list
  return(list(edge_list = edge_list, TOM = TOM))
}

# Function to write the edge list in a csv
writeEdgeList <- function(eset_file,netwk_type,edge_list){
  # Replace spaces with underscores in netwk_type
  networkType <- gsub(" ", "_", netwk_type)
  
  # Extract the dataset ID we are studying:
  dataset_id <- gsub("_eset.rds", "", eset_file)
  
  # Build the filename based on network type and dataset ID without spaces
  filename = sprintf("~/Internship_data/CSU/WGCNA/Results/edgelist_%s_%s.tsv",dataset_id,networkType)
  
  # Export the network file to be read into Cytoscape, VisANT, etc.
  write_delim(edge_list, 
              file = filename, 
              delim = "\t")
}

# Function to write edge list but for each module only not for all 
writeEdgeListForModule <- function(edge_list, module, netwk_type, eset_file) {
  # Replace spaces with underscores in netwk_type
  networkType <- gsub(" ", "_", netwk_type)
  # Extract the dataset ID we are studying:
  dataset_id <- gsub("_eset.rds", "", eset_file)
  
  # Construire le nom du fichier basé sur le module et le type de réseau
  filename <- sprintf("~/Internship_data/CSU/WGCNA/Results/edgelist_%s_%s_%s.tsv",dataset_id,networkType,module)
  
  # Exporter la liste des arêtes pour être lue dans Cytoscape, VisANT, etc.
  write_delim(edge_list, file = filename, delim = "\t")
}


# Function to create an interactive table
createTable <- function(df){
  display_table(df,
                interactive = TRUE,
                highlighted_color = "yellow",
                nmax_display = 10)
}

# Function to choose hub genes in each module
chooseTopHubGenes <- function(list_of_modules, exprs_data, power, netwk_type,omitColors=F) {
  top_hub_genes <- list()
  
  hub_genes_df <- data.frame(Module = character(), HubGene = character(), stringsAsFactors = FALSE)
  
  for (module in names(list_of_modules)) {
    module_genes <- list_of_modules[[module]]
    if (length(module_genes) > 1) {
      exprs_module <- exprs_data[,module_genes, drop = FALSE]
      top_hub_genes[[module]] <- chooseTopHubInEachModule(
        exprs_module,
        colorh = rep(module, ncol(exprs_module)),
        power = power,
        type = netwk_type,
        omitColors = omitColors
      )
    }
  }
  
  return(top_hub_genes)
}