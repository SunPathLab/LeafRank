#' calculate the up messages along each branch to the up node
#' branch index * up node state
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param time_scale time normalization factor
#' @param argument a list of birth, death and mutation data
#' @param E_list pre-calculated E_list (see the next step)
#' @param non_negativity_cutoff
#' @return matrix of "messages" sent to a node by its progeny lineages, rows: node; columns: fitness type
#' @importFrom Brobdingnag as.brobmat
#' @export
calc_up_messages <- function(phy, time_scale, argument, E_list, T_vector, non_negativity_cutoff){
  
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  
  ## check if the up message from a node has been calculated
  cal_status <- numeric(num_of_node)
  
  b <- argument[[1]]
  fitness_count <- length(b)
  
  up_messages <-  replicate(fitness_count, numeric(num_of_branch))
  
  time_estimate <- node_time_to_present(phy, time_scale)
  
  jobs_done <- 0
  while (jobs_done < num_of_branch) {
    print(jobs_done)
    for (i in 1:num_of_branch) {
      ## print(i)
      up_node <- edge_data[i,1]
      down_node <- edge_data[i,2]
      t_1 <- time_estimate[down_node]
      #print(t_1)
      
      if (cal_status[down_node]==0){
        
        cal_ready <- 1   # if the other required messages are calculated
        ## if it is internal
        if (node_status(phy, down_node)){
          for (j in 1:num_of_branch){
            if (edge_data[j,1]==down_node){
              cal_ready <- cal_ready * cal_status[edge_data[j,2]]
            }
          }
        }
        
        ## if the branch connects to an external node
        if (!node_status(phy, down_node)) {

          # scale the time to represent the actual time of growth of the tumor
          t <- length_data[i]/time_scale
          
          ######
          ## i: branch number
          ## each element represents a state of up node
          up_messages[i,] <- log(rowSums(integrate_prop(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff)))
          ######
          
          jobs_done <- jobs_done + 1
          cal_status[down_node] <- 1
        } else if (cal_ready) {   #internal, ready to be calculated
          
          
          t <- length_data[i]/time_scale
          
          ######
          ## each element represents a state of down node
          ## i: branch index
          ## down_node: index of down node
          temp_1 <- up_m_des(up_messages, down_node, phy, argument)
          ######
          
          ######
          ## up node * down node
          temp_2 <- integrate_prop(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff)
          ######
          
          temp_1 <- Brobdingnag::as.brobmat(temp_1)
          temp_2 <- Brobdingnag::as.brobmat(temp_2)

          ## temp_3 final
          for (k in 1:fitness_count){
            temp_3 <- 0
            for (l in 1:fitness_count){
              temp_3 <- temp_3 + temp_2[k,l]*temp_1[l]
            }
            up_messages[i,k] <- log(temp_3)
          }
          jobs_done <- jobs_done + 1
          cal_status[down_node] <- 1

        }
      }
    }
  }
  return(up_messages)
}


#' calculate the down messages 
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param time_scale time normalization factor
#' @param argument a list of birth, death and mutation data
#' @param E_list pre-calculated E_list (see the next step)
#' @param up_messages pre-calculated up messages
#' @param non_negativity_cutoff
#' @return matrix of "messages" sent to a node by its ancestral and sibling lineages, rows: node; columns: fitness type
#' @importFrom Brobdingnag as.brob
#' @export
calc_down_messages <- function(phy, time_scale, argument, E_list, up_messages, T_vector, non_negativity_cutoff){
  
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  
  ## check if the down message to a node has been calculated
  cal_status <- numeric(num_of_node)
  
  b <- argument[[1]]
  fitness_count <- length(b)
  
  ## num_of_branch by fitness_count zero matrix
  down_messages <-  replicate(fitness_count, numeric(num_of_branch))
  
  ## time to present estimate of each node
  time_estimate <- node_time_to_present(phy, time_scale)
  
  ## get the id number for the root
  root_id <- 0
  for (i in 1:num_of_node){
    if (!is.element(i, edge_data[,2])){
      root_id <- i
    }
  }
  
  jobs_done <- 0
  
  
  while (jobs_done < num_of_branch){
    ## print(jobs_done)
    for (i in 1:num_of_branch) {
      print(i)
      up_node <- edge_data[i,1]
      down_node <- edge_data[i,2]
      t_1 <- time_estimate[down_node]
      cal_ready <- 0
      
      if (cal_status[up_node]==1 | up_node == root_id) {
        cal_ready <- 1  ## ready if the down messages to the up node is ready or it is the root node
      }
      
      if (cal_status[down_node]==0 & cal_ready == 1) {
        
        
 
        t <- length_data[i]/time_scale
        temp_1 <- integrate_prop(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff)



        if (up_node!= root_id){
          temp_2 <- up_m_sib(up_messages, down_messages, up_node, phy, argument, i)
        } else {
          temp_2 <- numeric(fitness_count) 
          temp_2 <- Brobdingnag::as.brob(temp_2)
          
          sibling_data <- get_sibling(phy)
          sibling_branch <- sibling_data[i]
          
          temp_x <- exp(Brobdingnag::as.brob(up_messages[sibling_branch,1]))*2*b[1]
          temp_y <- exp(Brobdingnag::as.brob(up_messages[sibling_branch,2]))*mu[1]
          temp_2[1] <- temp_x + temp_y
          
          temp_y <- exp(Brobdingnag::as.brob(up_messages[sibling_branch,1]))*mu[1]
          temp_2[2] <- temp_y
        }

        for (k in 1:fitness_count) {

          # temp_3 final output
          temp_3 <- 0
          for (l in 1:fitness_count) {
            temp_3 <- temp_3 + temp_2[l]*temp_1[l,k]
          }
          down_messages[i,k] <- log(temp_3)
          ## print(down_messages[i,k])
        }
        
        ## down_messages[i,] <- log(temp_2 %*% temp_1)
        jobs_done <- jobs_done + 1
        cal_status[down_node] <- 1
      }
    }
  }             
  return(down_messages)
}


#' obtain the message upcoming from the down_node by taking into account its branching
#' each element represents a state of down node of the branch of interest
#' handle the branching from the child
#'
#' 
#' @param up_messages pre-calculated up messages
#' @param down_node index of the down node
#' @param phy a phylo object (tree object returned by nj)
#' @param argument a list of birth, death and mutation data
#' @return the message (across fitness states) upcoming from the down node
#' @importFrom Brobdingnag as.brob
#' @export
up_m_des <- function(up_messages, down_node, phy, argument) {
  b <- argument[[1]]
  mu <- argument[[3]]
  fitness_count <- length(b)
  
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  
  ## get the index of two descendants branch of the down_node (children of the child)
  des <- c()
  for (j in 1:num_of_branch) {
    if (edge_data[j,1]==down_node) {
      des <- c(des, j) 
    }
  }
  index_1 <- des[1]
  index_2 <- des[2]
  
  sol <- replicate(fitness_count,0)
  sol <- Brobdingnag::as.brob(sol)
  
  for (i in 1:fitness_count-1) {
    temp_1 <- up_messages[index_1,i]+up_messages[index_2,i]+log(2*b[i])
    temp_2 <- up_messages[index_1,(i+1)]+up_messages[index_2,i]+log(mu[i])
    temp_3 <- up_messages[index_1,i]+up_messages[index_2,(i+1)]+log(mu[i])
    temp_1 <- exp(Brobdingnag::as.brob(temp_1))
    temp_2 <- exp(Brobdingnag::as.brob(temp_2))
    temp_3 <- exp(Brobdingnag::as.brob(temp_3))
    sol[i] <- temp_1 + temp_2 + temp_3
  }

  
  i <- fitness_count
  temp_1 <- up_messages[index_1,i]+up_messages[index_2,i]+log(2*(b[i]+mu[i]))
  temp_1 <- exp(Brobdingnag::as.brob(temp_1))
  sol[i] <- temp_1
  
  
  return (sol)  
}


#' obtain the message upcoming from the up_node by taking into account the branching of its parental node
#' for calculating down message, we need to know the branching staring the mother node (sibling lineages)
#' convert numbers in log form to a probability measure
#' 
#' 
#' @param up_messages pre-calculated up messages
#' @param down_messages pre-calculated down messages
#' @param up_node index of the up node
#' @param phy a phylo object (tree object returned by nj)
#' @param argument a list of birth, death and mutation data
#' @param i branch index
#' @return the message (across fitness states) coming from the up_node
#' @importFrom Brobdingnag as.brob
#' @export
up_m_sib <- function(up_messages, down_messages, up_node, phy, argument, i){
  b <- argument[[1]]
  mu <- argument[[3]]
  fitness_count <- length(b)
  
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  
  ## get the index of parent branch
  for (j in 1:num_of_branch){
    if (edge_data[j,2]==up_node){
      parent_branch <- j
    }
  }

  ## index of the sibling branch
  sibling_data <- get_sibling(phy)
  sibling_branch <- sibling_data[i]
  
  sol <- replicate(fitness_count,0)
  sol <- Brobdingnag::as.brob(sol)
  
  if (fitness_count>2){
    for (i in 2:(fitness_count-1)) {
      temp_1 <- up_messages[sibling_branch,i]+down_messages[parent_branch,i]+log(2*b[i])
      temp_2 <- up_messages[sibling_branch,(i+1)]+down_messages[parent_branch,i]+log(mu[i])
      temp_3 <- up_messages[sibling_branch,(i-1)]+down_messages[parent_branch,(i-1)]+log(mu[i])
      temp_1 <- exp(Brobdingnag::as.brob(temp_1))
      temp_2 <- exp(Brobdingnag::as.brob(temp_2))
      temp_3 <- exp(Brobdingnag::as.brob(temp_3))
      sol[i] <- temp_1 + temp_2 + temp_3
    }
  }
  

  ## target branch starts at state 1
  i <- 1
  temp_1 <- up_messages[sibling_branch,i]+down_messages[parent_branch,i]+log(2*b[i])
  temp_1 <- exp(Brobdingnag::as.brob(temp_1))
  temp_2 <- up_messages[sibling_branch,(i+1)]+down_messages[parent_branch,i]+log(mu[i])
  temp_2 <- exp(Brobdingnag::as.brob(temp_2))
  sol[i] <- temp_1 + temp_2

  ## target branch starts at the highest state
  i <- fitness_count
  temp_1 <- up_messages[sibling_branch,i]+down_messages[parent_branch,i]+log(2*b[i])
  temp_1 <- exp(Brobdingnag::as.brob(temp_1))
  temp_2 <- up_messages[sibling_branch,(i-1)]+down_messages[parent_branch,(i-1)]+log(mu[i])
  temp_2 <- exp(Brobdingnag::as.brob(temp_2))
  sol[i] <- temp_1 + temp_2
  
  return (sol)  
}



#' convert numbers in log form to a probability measure
#' 
#' @param prob_vector a vector of probablity-related values that need to be normalized
#' @return a vector of normalized probablities
#' @export
log_normalize <- function(prob_vector) {
  prob_vector <- prob_vector - max(prob_vector)
  prob_vector <- exp(prob_vector)
  prob_vector <- prob_vector/sum(prob_vector)
  return(prob_vector)
}



#' calculate the marginal distribution for each node
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param time_scale time normalization factor
#' @param up_messages pre-calculated "up-messages" to every nodes
#' @param down_messages pre-calculated "down-messages" to every nodes
#' @param argument a list of birth, death and mutation data
#' @return a matrix of marginal probablity distribution on fitness type (column) for each node (row)
#' @importFrom Brobdingnag as.brob
#' @export 
calc_marginal_probabilities <- function(phy, time_scale, up_messages, down_messages, argument){
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  b <- argument[[1]]
  fitness_count <- length(b)
  mu <- argument[[3]]
  marginal_prob <-  replicate(fitness_count, numeric(num_of_node))
  
  root_id <- 0
  for (i in 1:num_of_node){
    if (!is.element(i, edge_data[,2])){
      root_id <- i
    }
  }
  
  for (i in 1:num_of_node){
    
    if (i != root_id){  ## if i is not the root node
      
      temp_1 <- numeric(fitness_count)
      for (j in 1:num_of_branch){
        if (edge_data[j,2]==i){
          temp_1 <- temp_1 + down_messages[j,]  ## get down messages to the target node
        }
      }
      
      branch_index <- c()  ## get two descendant branches for the target node (empty set for external nodes)
      for (j in 1:num_of_branch){
        if (edge_data[j,1]==i){
          branch_index <- c(branch_index, j)
        }
      }
      branch_1 <- branch_index[1]
      branch_2 <- branch_index[2]
      
      temp_2 <- numeric(fitness_count)
      temp_2 <- Brobdingnag::as.brob(temp_2)
      if (!is.null(branch_1)){  ## internal nodes
        
        for (k in 1:fitness_count){
          if (k==fitness_count){
            temp_2[k] <- exp(Brobdingnag::as.brob(up_messages[branch_1, k])) * exp(Brobdingnag::as.brob(up_messages[branch_2, k])) * 2 * (b[k]+mu[k])
          } else {
            temp_2[k] <- exp(Brobdingnag::as.brob(up_messages[branch_1, k])) * exp(Brobdingnag::as.brob(up_messages[branch_2, k])) * 2 * b[k]
            temp_2[k] <- temp_2[k] + exp(Brobdingnag::as.brob(up_messages[branch_1, (k+1)])) * exp(Brobdingnag::as.brob(up_messages[branch_2, k])) * mu[k]
            temp_2[k] <- temp_2[k] + exp(Brobdingnag::as.brob(up_messages[branch_1, k])) * exp(Brobdingnag::as.brob(up_messages[branch_2, (k+1)])) * mu[k]
          }  ## get up messages to the target node
        }
        temp_2 <- log(temp_2)
        marginal_prob[i,] <- log_normalize(temp_1+temp_2)
        
      } else {

        marginal_prob[i,] <- log_normalize(temp_1)  ## for external nodes, we only have down messages
      }
    } else {  ## root node, this part can be simplified

      temp_1 <- numeric(fitness_count)
      temp_1[1] <- 1
      
      branch_index <- c()
      for (j in 1:num_of_branch){
        if (edge_data[j,1]==i){
          branch_index <- c(branch_index, j)
        }
      }
      branch_1 <- branch_index[1]
      branch_2 <- branch_index[2]
      
      temp_2 <- numeric(fitness_count)
      temp_2 <- Brobdingnag::as.brob(temp_2)
      
      j <- 1

      temp_2[j] <- exp(Brobdingnag::as.brob(up_messages[branch_1, j])) * exp(Brobdingnag::as.brob(up_messages[branch_2, j])) * 2*b[j]
      temp_2[j] <- temp_2[j] + exp(Brobdingnag::as.brob(up_messages[branch_1, (j+1)])) * exp(Brobdingnag::as.brob(up_messages[branch_2, j])) * mu[j]
      temp_2[j] <- temp_2[j] + exp(Brobdingnag::as.brob(up_messages[branch_1, j])) * exp(Brobdingnag::as.brob(up_messages[branch_2, (j+1)])) * mu[j]
 
      temp_2 <- log(temp_2)
      marginal_prob[i,] <- log_normalize(temp_1+temp_2)
      }
      
      
      
    }
  
  return(marginal_prob)
}


#' calculate the marginal distribution for each node
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param marginal_prob a matrix of marginal probablity distribution on fitness type (column) for each node (row)
#' @param argument a list of birth, death and mutation data
#' @return vector of mean fitness of each node  
#' @export 
mean_fitness <- function(phy, marginal_prob, argument){
  edge_data <- phy$edge
  length_data <- phy$edge.length
  num_of_node <- node_num(phy)
  num_of_branch <- dim(edge_data)[1]
  b <- argument[[1]]
  fitness_count <- length(b)
  
  mean_result <- numeric(num_of_node)
  for (i in 1:num_of_node){
    for (j in 1:fitness_count){
      mean_result[i] <- mean_result[i] + b[j]*marginal_prob[i,j]
    }
  }
  return(mean_result)
}

