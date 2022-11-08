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
    phy <- ape::nj(M)
    phy.rooted <- ape::root(phy, outgroup = "normal", resolve.root = TRUE)
    phy <- ape::drop.tip(phy.rooted, "normal", trim.internal = TRUE, subtree = FALSE)
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
#' @param mu driver mutation rates for types (vector)
#' @param T_vector time discretization (vector)
#' @param non_negativity_cutoff 
#' @return the mean fitness vector
#' @export
ith.Fitness <- function(phy, outFile, rho, d_t, time_scale, b_rate, d_rate, mu, T_vector, non_negativity_cutoff){

    argument = list(b_rate, d_rate, mu)
    
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

    # up messages
    message("calculating up messages")
    up_messages=calc_up_messages(phy, time_scale, argument, E_list, T_vector, non_negativity_cutoff)

    # down messages
    message("calculating down messages")
    down_messages=calc_down_messages(phy, time_scale, argument, E_list, up_messages, T_vector, non_negativity_cutoff)

    # marginal probabilities
    message("calculating marginal probabilities")
    marginal_prob <- calc_marginal_probabilities(phy, time_scale, up_messages, down_messages, argument)

    # mean fitness
    message("calculating mean fitness")
    mean_result <- mean_fitness(phy, marginal_prob, argument)

    # saveResult
    saveRDS(mean_result, file=outFile)

    return("success!")
}

