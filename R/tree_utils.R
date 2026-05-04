#' return the number of nodes (internal and external)
#'
#' @param phy a phylo object (tree object returned by nj)
#' @return the number of nodes in the tree
#' @export
node_num <- function(phy){
  num_of_node <- max(phy$edge)
  return (num_of_node)  
}


#' Check if a node is internal or external
#'
#' @param phy a phylo object (tree object returned by nj)
#' @param node_id the id of the node of interest
#' @return (logical) TRUE: internal, FALSE: external
#' @export
node_status <- function(phy, node_id){
  edge_data <- phy$edge
  result <- is.element(node_id, edge_data[,1])
  return (result)  
}



#' Return the depth of each nodes (depth measures the distance between the root and a node)
#'
#' @param phy a phylo object (tree object returned by nj)
#' @return vector of depths of nodes
#' @export
get_depth <- function(phy){
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  depth_data <- numeric(num_of_node)
  depth_done <- numeric(num_of_node)
  
  ## get the id number for the root
  root_id <- 0
  for (i in 1:num_of_node){
    if (!is.element(i, edge_data[,2])){
      root_id <- i
    }
  }
  
  depth_done[root_id] <- 1
  jobs_done <- 1
  while (jobs_done < num_of_node){
    for (i in 1:num_of_branch){
      up_node <- edge_data[i,1]
      down_node <- edge_data[i,2]
      if (depth_done[up_node]==1 & depth_done[down_node]==0){
        jobs_done <- jobs_done + 1
        depth_data[down_node] <- depth_data[up_node] + length_data[i]
        depth_done[down_node] <- 1
      }
    }
  }
  return(depth_data)
}



#' get the ids of external nodes that are descendants of an internal node
#' return null if the input node is an external node
#'
#' @param phy a phylo object (tree object returned by nj)
#' @param node_id the id of the node of interest
#' @return vector ids of the two children nodes
#' @export
get_descendant <- function(phy, node_id){
  edge_data <- phy$edge
  num_of_branch <- dim(edge_data)[1]
  descendants <- c()
  if (node_status(phy, node_id)){
    parents <- c(node_id)
    while (length(parents)>0){
      remove_list <- c()
      add_list <- c()
      for (i in 1:num_of_branch){
        up_node <- edge_data[i,1]
        down_node <- edge_data[i,2]
        if (is.element(up_node, parents) & (!node_status(phy, down_node))){
          descendants <- c(descendants, down_node)
          remove_list <- c(remove_list, up_node)
        } else if (is.element(up_node, parents) & (node_status(phy, down_node))){
          add_list <- c(add_list, down_node)
          remove_list <- c(remove_list, up_node)
        }
      }
      parents = parents[!(parents %in% remove_list)]
      parents <- c(parents, add_list)
    }
    return (descendants)
  } else{
    return (descendants) 
  }
}



#' transform the tree into an ultrametric tree
#' re-scale the branch length (measured by mutations) to the real time scale
#' 
#' time = 0 for leaves
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param time_scale time normalization factor
#' @return vector of normalized branch lengths in real time scale
#' @export
node_time_to_present <- function(phy, time_scale){
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  time_estimate <- numeric(num_of_node)
  depth_data <- get_depth(phy) 
  for (i in 1:num_of_node){
    if (!node_status(phy, i)) {
      time_estimate[i] <- 0
    } else {
      descendants <- get_descendant(phy, i)
      for (j in descendants){
        time_estimate[i] <- time_estimate[i] + depth_data[j]
      }
      time_estimate[i] <- (time_estimate[i]/length(descendants)-depth_data[i])/time_scale
    }
  }
  return (time_estimate)  
}


#' @param phy a phylo object (tree object returned by nj)
#' @return vector of branch ids corresponding to the sibling branch (of ordered branch ids)
#' @export
get_sibling <- function(phy){
  edge_data <- phy$edge
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  sibling_data <- numeric(num_of_branch)
  for (i in 1:num_of_branch){      # for each branch id
    up_node <- edge_data[i,1]
    down_node <- edge_data[i,2]
    for (j in 1:num_of_branch) {
      if ((edge_data[j,1]==up_node)&(edge_data[j,2]!=down_node)){   # get the sibling branch id
        sibling_data[i] <- j
      }
    }
  }
  return (sibling_data)  
}


#' @param tree a phylo object with feature node_wgd, indicating the WGD status of tip nodes.
#' @return a phylo object with internal nodes and branches labeled accordingly for WGD status.
#' @export
label_internal_nodes <- function(tree) { 
  num_branch <- length(tree$edge[,1])
  tree$edge_wgd <- rep(-1, num_branch)
  job_done <- 0
  while (job_done < num_branch){
    for (i in 1:num_branch){
      if (tree$edge_wgd[i] != -1 ){
        next
      }
      child <- tree$edge[i,2]
      parent <- tree$edge[i,1]
      if (tree$node_wgd[child] == -1){
        next
      }
      if (tree$node_wgd[child] == 1){
        tree$edge_wgd[i] = 1
        job_done = job_done + 1
      }else{
        tree$edge_wgd[i] = 0
        job_done = job_done + 1
      }
      if (tree$node_wgd[parent] == -1){
        sib <- get_siblings(child,tree)
        if (tree$node_wgd[sib] == -1){
          next
        }else if(tree$node_wgd[sib] == tree$node_wgd[child]){
          tree$node_wgd[parent] = tree$node_wgd[child]
        }else{
          tree$node_wgd[parent] = 0
        }
      }
      
    }
  }
  tree
}

#' @param tree a phylo object with feature node_wgd, indicating the WGD status of tip nodes.
#' @return a phylo object with internal nodes and branches labeled accordingly for WGD status based on the Des. assumption.
#' @export
label_internal_nodes_des <- function(tree) { 
  num_branch <- length(tree$edge[,1])
  tree$edge_wgd <- rep(-1, num_branch)
  job_done <- 0
  while (job_done < num_branch){
    for (i in 1:num_branch){
      if (tree$edge_wgd[i] != -1 ){
        next
      }
      child <- tree$edge[i,2]
      parent <- tree$edge[i,1]
      if (tree$node_wgd[child] == -1){
        next
      }
      if (tree$node_wgd[parent] == 1){
        tree$edge_wgd[i] = 1
        job_done = job_done + 1
      }else{
        tree$edge_wgd[i] = 0
        job_done = job_done + 1
      }
      if (tree$node_wgd[parent] == -1){
        sib <- get_siblings(child,tree)
        if (tree$node_wgd[sib] == -1){
          next
        }else if(tree$node_wgd[sib] == tree$node_wgd[child]){
          tree$node_wgd[parent] = tree$node_wgd[child]
        }else{
          tree$node_wgd[parent] = 0
        }
      }
      
    }
  }
  tree
}

#' Get a full set of LeafRank configuration. 
#'
#' @param phy a phylo object (tree object returned by nj)
#' @return vector of branch ids corresponding to the sibling branch (of ordered branch ids)
#' @export
get_full_pars <- function(
    b_rates = NULL,
    d_rates = NULL,
    nu      = NULL,
    rho     = NULL,
    tau     = NULL,
    init    = NULL,
    tree,
    model   = 'default'
  ){
  cell_num <- length(tree$tip.label)
  input <- list(b_rates, d_rates, nu, rho, tau)
  input_A <- list(b_rates,d_rates)
  if (all(sapply(input_A, is.null))){
    num_pheno <- 16
  }else if(is.null(input_A)){
    idx <- which(!sapply(input_A,is.null))
    num_pheno <- length(input_A[[idx]])
  }else{
    num_pheno <- length(b_rates)
  }
  p_init <- matrix(rep(0,num_pheno),nrow = 1)
  p_init[1] <- 1
  if (is.null(init)){
    init <- p_init
  }else if(length(init) < num_pheno){
    init <- p_init
  }
  sum_vec <- matrix(rep(1,num_pheno),ncol = 1)
  
  b_coeff <- 4^(1/(num_pheno-1))
  p_b_rates <- 1.1*b_coeff^(0:num_pheno-1)-1.1
  p_d_rates <- replicate(num_pheno, 1)
  p_nu      <- 0.0001
  p_rho     <- 0.0005
  
  while (sum(sapply(input, is.null))> 1) {
    idx <- which(sapply(input, is.null))
    if (idx[1] == 1){
      if (idx[2] >2){
        b_rates <- d_rates + p_b_rates
      }else{
        b_rates <- 1.1 + p_b_rates
      }
    }else if(idx[1] == 2){
      d_rates <- p_d_rates
    }else if(idx[1] == 3){
      nu <- p_nu
    }else if(idx[1] == 4){
      rho <- p_rho
    }
    input <- list(b_rates, d_rates, nu, rho, tau)
  }
  
  idx <- which(sapply(input, is.null))
  
  if (idx == 5){
    total_cell <- cell_num/rho
    A <- diag(b_rates-d_rates) + cbind(matrix(0,ncol = 1,nrow = num_pheno),diag(nu,nrow = num_pheno, ncol = num_pheno-1))
    f <- function(t){
      as.numeric(init %*% expm::expm(A*t) %*% sum_vec)-total_cell
    }
    z_sign <- sign(f(0))
    t <- 0
    for (i in c(10,20,50,100,200,500)){
      if (sign(f(i))!= z_sign){
        res <- uniroot(f,interval  = c(0,i))
        t <- res$root
        break
      }
    }
    if (t == 0){
      stop("Expected time is out of default range: Need manurally extend the range")
    }
    tau <- 1/t
  }else if (idx == 4){
    A <- diag(b_rates-d_rates) + cbind(matrix(0,ncol = 1,nrow = num_pheno),diag(nu,nrow = num_pheno, ncol = num_pheno-1))
    total_cell <- as.numeric(init %*% expm::expm(A*(1/tau)) %*% sum_vec)
    rho <- cell_num/total_cell
  }
  
  input <- list(b_rates, d_rates, nu, rho, tau)
  return(input)
}

