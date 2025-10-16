

# set.seed(1)
# outputMCL <- MakeClusterList(source, 1, 1)
# save(outputMCL[[1]], file="ex_common_clusters.rda")

  #CHANGE ME IF NEEDED
  ptmtable <- ex_tiny_ptm_table 

  set.seed(1)
  OutputMCL <- MakeClusterList(ptmtable)
  
  ex_common_cluster         <- OutputMCL[[1]]
  ex_adj_consensus          <- OutputMCL[[2]]
  ex_ptm_correlation_matrix <- OutputMCL[[3]]
  #save
  
  
  
  OutputMCN <- MakeCorrelationNetwork(ex_adj_consensus, ex_ptm_correlation_matrix)
  ex_ptm_cccn_g      <- OutputMCN[[1]]
  ex_gene_cccn_g     <- OutputMCN[[2]]
  ex_ptm_cccn_edges  <- OutputMCN[[3]]
  ex_gene_cccn_edges <- OutputMCN[[4]]
  #save

  
  #TO DO import bioplanet file
  OutputPCN <- PathwayCrosstalkNetwork(ex_common_clusters, bioplanet.file = "pathway.csv", createfile = FALSE)
  ex_pathway_crosstalk_network <- PCN.data[[1]]
  ex_PCNedgelist <- PCN.data[[2]]
  ex_pathways_list <- PCN.data[[3]]
  #save
  
  
  stringdb.edges <- GetSTRINGdb(ex_gene_cccn_edges)
  MakeDBInput(ex_gene_cccn_edges, file.path.name = "db_nodes.txt")
  kinsub.edges <- format.kinsub.table(kinasesubstrate.filename = "Kinase_Substrate_Dataset.txt")
  
  
  
  
  