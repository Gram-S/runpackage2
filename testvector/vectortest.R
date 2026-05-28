  LowSolution <- function(clusters.list, keeplength=3) { # >>>> NEW method
    
    mat <- matrix(0, nrow = length(PTMnames), ncol = length(PTMnames), dimnames = list(PTMnames, PTMnames))
    adjacency_matrices <- list(mat, mat, mat)
    count <- 0
    
    for (cluster_group in clusters.list){
      count <- count + 1
      sapply(cluster_group, function(x) adjacency_matrices[[count]][ x[[1]], x[[1]] ] <<- 1 )
    }
    
    return(adjacency_matrices)
  }
  
  CurrentSolution <- function(clusters.list, keeplength=3) { # >>>> NEW method
    
    co_membership_matrix <- function(clusters, all_ptms) {
      # 1. square matrix of all PTMs
      mat <- matrix(0, nrow = length(PTMnames), ncol = length(PTMnames),
                    dimnames = list(PTMnames, PTMnames))
      # 2. For each cluster, set all PTM–PTM pairs in the cluster to 1 (indicating co-membership).
      for (cluster in clusters) {
        ptms <- cluster$PTMnames
        if (length(ptms) > 1) {
          mat[ptms, ptms] <- 1
        }
      }
      diag(mat) <- 0
      return(mat)
    }
    
    adjacency_matrices <- lapply(clusters.list, co_membership_matrix)
  
    return(adjacency_matrices)
  }
  
  ReduceSolution <- function(clusters.list, keeplength=3) { # >>>> NEW method
    
    co_membership_matrix <- function(clusters, all_ptms) {
      # 1. square matrix of all PTMs
      mat <- matrix(0, nrow = length(PTMnames), ncol = length(PTMnames),
                    dimnames = list(PTMnames, PTMnames))
      # 2. For each cluster, set all PTM–PTM pairs in the cluster to 1 (indicating co-membership).
      for (cluster in clusters) {
        ptms <- cluster$PTMnames
        if (length(ptms) > 1) {
          mat[ptms, ptms] <- 1
        }
      }
      diag(mat) <- 0
      return(mat)
    }
    
    adjacency_matrices <- purrr::map(clusters.list, co_membership_matrix, all_ptms=all_ptms)
    message("REDUCE SOLUTION for common clusters takes: ", Sys.time() - start_time)
    
    return(adjacency_matrices)
  }
  
  low_time <- Sys.time()
  l <- LowSolution(clusters.list)
  message("LOW SOLUTION for common clusters takes: ", Sys.time() - low_time)
  
  current_time <- Sys.time()
  c <- CurrentSolution(clusters.list) 
  message("CURRENT SOLUTION for common clusters takes: ", Sys.time() - current_time)
 
  reduce_time <- Sys.time()
  r <- ReduceSolution(clusters.list)
  message("REDUCE SOLUTION for common clusters takes: ", Sys.time() - reduce_time)
  
  print(all.equal(c, r))
  

