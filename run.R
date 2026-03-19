# How to run
# htmlwidgets::saveWidget(profvis::profvis({codehere}), file="Profiling/name") # DO NOT APPEND .html TO FILENAME 
# ---



#CHANGE ME IF NEEDED
  ex_ptmtable <- ex_tiny_ptm_table 
  ex_bioplanet <- ex.bioplanet

  set.seed(88)
  
  #MakeClusterList
  OutputMCL <- MakeClusterList(ex_ptmtable)
  ex_common_clusters        <- OutputMCL[[1]]
  ex_adj_consensus          <- OutputMCL[[2]]
  ex_ptm_correlation_matrix <- OutputMCL[[3]]

  
  
  #MakeCorrelationNetwork
  OutputMCN <- MakeCorrelationNetwork(ex_adj_consensus, ex_ptm_correlation_matrix)
  ex_ptm_cccn_edgelist  <- OutputMCN[[1]]
  ex_gene_cccn_edgelist <- OutputMCN[[2]]
  ex_genelist           <- OutputMCN[[3]]
  
  
  
  #Get stringdb edges
  ex_stringdb_edges <- GetSTRINGdb.edges(ex_gene_cccn_edgelist, ex_genelist)
  
  #Make Genemania Input 
  ex_genemania_edges <- MakeDBInput(ex_gene_cccn_edgelist, file.path.name = "db_nodes.txt")
  
  #Get Kinsub edges - TO DO 
  #ex_kinsub_edges <- formatKinsubTable(kinasesubstrate.filename = "phospho_cleaned_mapped.txt", ex_gene_cccn_edgelist) #WHY ARE THE PARAMETERS LIKE THIS WHYYYYYY
  
  
  
  #BuildClusterFilteredNetwork
  OutputBCFN <- BuildClusterFilteredNetwork(ex_gene_cccn_edgelist, ex_stringdb_edges, ex_genemania_edges, NULL)
  ex_combined_ppi <- OutputBCFN[[1]]
  ex_cfn          <- OutputBCFN[[2]]


  #BuildPathwayCrosstalkNetwork
  OutputBPCN <- BuildPathwayCrosstalkNetwork(ex_common_clusters, ex_bioplanet, createfile = FALSE)
  ex_pathway_crosstalk_network <- OutputBPCN[[1]]
  ex_PCNedgelist <- OutputBPCN[[2]]
  ex_pathways_list <- OutputBPCN[[3]]

  
  # Save all the example data
   # for(name in ls(.GlobalEnv)){
   #  save(name, file=paste(name, '.rda', sep=''))
   # }
  
  
  
