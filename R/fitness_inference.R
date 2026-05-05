#‘ Apply the neighbor joining method to reconstruct the genealogical tree (phy)
#' prepare input of tree from mutation matrix
#' mutation input, prepare phy
#' the treeFile input (tab delimited) currently has this format
#' mutationID   cell1    cell2   cell3   cell4  cell5  ...
#' 12345          1        1       0       0      0    ...
#' 67890          0        1       1       1      0    ...
#' ...
#' 
#' @param treeFile the mutationMatrix file input (tab delimited)
#' @return a phylo object of the tree, with the root edge being removed
#' @importFrom ape nj 
#' @importFrom ape drop.tip
#' @importFrom ape root
#' @export
mutationMatrix2Tree <- function(treeFile) {
    treeData = read.delim(treeFile, stringsAsFactors = FALSE)
    samples = colnames(treeData)[-1]
    samples = c("normal", samples)
    treeData = data.frame(treeData, normal = 0)
    M <- matrix(0, length(samples), length(samples))
    dimnames(M) <- list(samples, samples)
    pairs = combn(length(samples), 2)
    for (c in 1:dim(pairs)[2]){
        sn1 = samples[pairs[1,c]]
        sn2 = samples[pairs[2,c]]
        M[sn1,sn2] = sum(treeData[,sn1] != treeData[,sn2])
        M[sn2,sn1] = sum(treeData[,sn2] != treeData[,sn1])
    }
    print(dim(M))
    phy <- ape::nj(M)
    phy.rooted <- ape::root(phy, outgroup = "normal", resolve.root = TRUE)
    phy <- ape::drop.tip(phy.rooted, "normal", trim.internal = TRUE, subtree = FALSE)
    return(phy)
}




#' prepare input of tree from distance matrix
#' distance matrix input, prepare phy
#' the distanceMatrix input (tab delimited) currently has this format
#'             cell1    cell2   cell3   cell4  cell5  ...
#' cell1          0        123     11       450     21    ...
#' cell2          123      0       34       34      45    ...
#' ...
#' the last cell represents the normal cell
#' @param distanceMatrix the mutationMatrix file input (comma delimited)
#' @return a phylo object of the tree, with the root edge being removed
#' @importFrom ape nj 
#' @importFrom ape drop.tip
#' @importFrom ape root
#' @export
distanceMatrix2Tree <- function(distanceMatrix) {
  M <- read.delim(distanceMatrix, header=F, sep=",")
  M <- as.numeric(unlist(M))
  matrixSize <- sqrt(length(M))
  M <- matrix(M, ncol = matrixSize, nrow = matrixSize)
  dimnames(M) <- list(1:matrixSize, 1:matrixSize)
  phy <- ape::nj(M)
  phy.rooted <- ape::root(phy, outgroup = matrixSize, resolve.root = TRUE)
  phy <- ape::drop.tip(phy.rooted, matrixSize, trim.internal = TRUE, subtree = FALSE)
  return(phy)
}



#' pipeline to do fitness inference
#' 
#' @param phy a phylo object (tree object returned by nj)
#' @param outFile output File with full path
#' @param rho sampling probability
#' @param d_t step size for integration
#' @param time_scale time normalization factor
#' @param b_rate birth rates for types (vector)
#' @param d_rate death rates for types (vector)
#' @param nu driver mutation rates for types (vector)
#' @param T_vector time discretization (vector)
#' @param non_negativity_cutoff Set threshold
#' @param use_parallel Indicator for parallel computing
#' @return the mean fitness vector
#' @export
LeafRank <- function(phy, outFile, rho, d_t, time_scale, b_rates, d_rates, nu, T_vector, non_negativity_cutoff, use_parallel = FALSE){
  
    argument = list(b_rates, d_rates, nu)

    #start calculation
    message("start: E_list")
    # E(t)
    E_sol <- integrate_phi_E(rho, T_vector, argument, d_t, non_negativity_cutoff)
    ## row: time, col: fitness
    E_list <- list()
    for (i in 1:dim(E_sol)[2]) {
        E_approx <- approxfun(T_vector, E_sol[,i], method="linear")    # approximate the functions of E_i(t) by interpolation
        E_list[[i]] <- E_approx
    }
    
    ## Compute the population integration
    edge_data <- phy$edge
    length_data <- phy$edge.length
    num_of_node <- node_num(phy)
    num_of_branch <- dim(edge_data)[1]
    
    time_estimate <- node_time_to_present(phy, time_scale) ## Approximated time of each node, with time scaled
    
    if (use_parallel){
      int_list <- foreach::foreach (i =1:num_of_branch, .packages = "LeafRank") %dopar% {
        up_node <- edge_data[i,1]
        down_node <- edge_data[i,2]
        t_1 <- time_estimate[down_node]
        t <- length_data[i]/time_scale
        temp <- integrate_prop(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff)
        temp
      }
    }else{
      message("Parallel packages not available, running sequentially.")
      int_list <- vector("list", num_of_branch)
      for (i in seq_len(num_of_branch)) {
        up_node <- edge_data[i,1]
        down_node <- edge_data[i,2]
        t_1 <- time_estimate[down_node]
        t <- length_data[i]/time_scale
        int_list[[i]] <- integrate_prop(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff)
      }
    }
    # up messages
    message("calculating up messages")
    up_messages=calc_up_messages(phy, time_scale, argument, rho, d_t, E_list, T_vector, non_negativity_cutoff, time_estimate, int_list)

    # down messages
    message("calculating down messages")
    down_messages=calc_down_messages(phy, time_scale, argument, rho, d_t, E_list, up_messages, T_vector, non_negativity_cutoff, time_estimate, int_list)

    # marginal probabilities
    message("calculating marginal probabilities")
    marginal_prob <- calc_marginal_probabilities(phy, time_scale, up_messages, down_messages, argument)

    # mean fitness
    message("calculating mean fitness")
    mean_result <- mean_fitness(phy, marginal_prob, argument)

    # saveResult
    outRes = list(phylo = phy, meanFitness = mean_result, upMessages = up_messages, downMessages = down_messages, marginalProb = marginal_prob)
    saveRDS(outRes, file=outFile)

    return("success!")
}





#' Compare the ranking inferred by two different sets of parameters
#'
#' @param result_1 mean fitness inferred by the first set of parameters
#' @param result_2 mean fitness inferred by the second set of parameters
#' @return a Kendall-tau like distance
#' @export
rank_diff <- function(result_1, result_2)
{
  x_1 <- length(result_1)
  x_1 <- (x_1+1)/2  # only check mean fitness for leaves
  kendall_dis <- 0
  
  for (i in 1:x_1){
    for (j in 1:x_1){
      if ((result_1[i]>result_1[j]) & (result_2[i]<result_2[j])) {
        kendall_dis <- kendall_dis + 1
        print(i)
        print(j)
      }
      if ((result_1[i]<result_1[j]) & (result_2[i]>result_2[j])) {
        kendall_dis <- kendall_dis + 1
        print(i)
        print(j)
      } 
    }
  }
  return (kendall_dis)
}



