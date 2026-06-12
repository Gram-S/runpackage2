# How to run
# htmlwidgets::saveWidget(profvis::profvis({codehere}), file="Profiling/name") # DO NOT APPEND .html TO FILENAME 
# ---

  # Relies on example data, so check if packaged is installed then load
  if(!require('PTMsToPathways')){
    stop('PTMsToPathways not found, please use devtools::document before running this')
  }
  
  # Load data from package
  data( 'ex_small_ptm_table', 'ex_bioplanet', package='PTMsToPathways' )

  #CHANGE ME IF NEEDED
  save_to_file <- FALSE # IF True then we write everything to the working directory
  ex_ptmtable <- ex_small_ptm_table 	
  set.seed(88)

  #MakeClusterList
  OutputMCL <- MakeClusterList(ex_ptmtable)
  ex_common_clusters        <- OutputMCL[[1]]
  ex_adj_consensus          <- OutputMCL[[2]]
  ex_ptm_correlation_matrix <- OutputMCL[[3]]

  
  
  #MakeCorrelationNetwork
  OutputMCN <- MakeCorrelationNetwork(ex_adj_consensus, ex_ptm_correlation_matrix)
  ex_ptm_cccn_edges  <- OutputMCN[[1]]
  ex_gene_cccn_edges <- OutputMCN[[2]]
  ex_gene_cccn_nodes <- OutputMCN[[3]]
  
  
  
  #Get stringdb edges
  ex_stringdb_edges <- GetSTRINGdb.edges(ex_gene_cccn_edges, ex_gene_cccn_nodes)
  
  #Make Genemania Input 
  # MakeDBInput(ex_gene_cccn_nodes, file.path.name = "ex_db_nodes.txt")
  ex_genemania_edges <- GetGeneMANIA.edges('ex_genemania_interactions.txt', ex_gene_cccn_nodes)
  
  #Get Kinsub edges - TO DO 
  #ex_kinsub_edges <- formatKinsubTable(kinasesubstrate.filename = "phospho_cleaned_mapped.txt", ex_gene_cccn_edgelist) #WHY ARE THE PARAMETERS LIKE THIS WHYYYYYY

  
  
  #BuildClusterFilteredNetwork
  OutputBCFN <- BuildClusterFilteredNetwork(ex_gene_cccn_edges, ex_stringdb_edges, ex_genemania_edges, NULL)
  ex_combined_ppi <- OutputBCFN[[1]]
  ex_cfn          <- OutputBCFN[[2]]


  #BuildPathwayCrosstalkNetwork
  OutputBPCN <- BuildPathwayCrosstalkNetwork(ex_common_clusters, ex_bioplanet, createfile = FALSE)
  ex_pathway_crosstalk_network <- OutputBPCN[[1]]
  ex_PCNedgelist <- OutputBPCN[[2]]
  ex_pathways_list <- OutputBPCN[[3]]

  rm('OutputMCL', 'OutputMCN', 'OutputBCFN', 'OutputBPCN')

  if(save_to_file){
  
  save(ex_common_clusters, file='ex_common_clusters.rda')
  save(ex_adj_consensus, file='ex_adj_consensus.rda')
  save(ex_ptm_correlation_matrix, file='ex_ptm_correlation_matrix.rda')

  save(ex_ptm_cccn_edges, file='ex_ptm_cccn_edges.rda')
  save(ex_gene_cccn_edges, file='ex_gene_cccn_edges.rda')
  save(ex_gene_cccn_nodes, file='ex_gene_cccn_nodes.rda')

  save(ex_stringdb_edges, file='ex_stringdb_edges.rda')
  save(ex_genemania_edges, file='ex_genemania_edges.rda')

  save(ex_combined_ppi, file='ex_combined_ppi.rda')
  save(ex_cfn, file='ex_cfn.rda')

  save(ex_pathway_crosstalk_network, file='ex_pathway_crosstalk_network.rda')
  save(ex_PCNedgelist, file='ex_PCNedgelist.rda')
  save(ex_pathways_list, file='ex_pathways_list.rda')

  }
  
  
  
