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
