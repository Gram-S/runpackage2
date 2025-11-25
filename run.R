  #CHANGE ME IF NEEDED
  ptmtable <- ex_tiny_ptm_table 

  set.seed(1)
  
  #MakeClusterList
  OutputMCL <- MakeClusterList(ptmtable)
  ex_common_clusters         <- OutputMCL[[1]]
  ex_adj_consensus          <- OutputMCL[[2]]
  ex_ptm_correlation_matrix <- OutputMCL[[3]]

  
  
  #MakeCorrelationNetwork
  OutputMCN <- MakeCorrelationNetwork(ex_adj_consensus, ex_ptm_correlation_matrix)
  ex_ptm_cccn_edgelist  <- OutputMCN[[1]]
  ex_gene_cccn_edgelist <- OutputMCN[[2]]
  ex_genelist           <- OutputMCN[[3]]
  
  
  
  #Get stringdb edges
  ex_stringdb_edges <- GetSTRINGdb(ex_gene_cccn_edgelist, ex_genelist)
  
  #Make Genemania Input 
  ex_genemania_edges <- MakeDBInput(ex_gene_cccn_edgelist, file.path.name = "db_nodes.txt")
  
  #Get Kinsub edges - TO DO 
  #ex_kinsub_edges <- formatKinsubTable(kinasesubstrate.filename = "phospho_cleaned_mapped.txt", ex_gene_cccn_edgelist) #WHY ARE THE PARAMETERS LIKE THIS WHYYYYYY
  
  
  
  #BuildClusterFilteredNetwork
  OutputBCFN <- BuildClusterFilteredNetwork(ex_gene_cccn_edgelist, ex_stringdb_edges, ex_genemania_edges, NULL)
  ex_combined_ppi <- OutputBCFN[[1]]
  ex_cfn          <- OutputBCFN[[2]]


  #PathwayCrosstalkNetwork
  OutputPCN <- PathwayCrosstalkNetwork(ex_common_clusters, bioplanet.file = "pathway.csv", createfile = FALSE)
  ex_pathway_crosstalk_network <- OutputPCN[[1]]
  ex_PCNedgelist <- OutputPCN[[2]]
  ex_pathways_list <- OutputPCN[[3]]
  
  
  
  
  #MakeDBInput(ex_gene_cccn_edges, file.path.name = "db_nodes.txt")
  #kinsub.edges <- format.kinsub.table(kinasesubstrate.filename = "Kinase_Substrate_Dataset.txt")
  
  
  
  
  