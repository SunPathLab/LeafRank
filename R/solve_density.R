#' Calculate Extinction density over time
#' note that E(t) is independent of tree structure
#' solve differential equations using Runge-Kutta methods (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
#'
#' @param f_initial vector of initial values (for each type of cells)
#' @param T_vector time discretization
#' @param argument a list of birth, death and mutation data (b <- argument[[1]]; d <- argument[[2]]; mu <- argument[[3]])
#' @param d_t step size for integration
#' @param non_negativity_cutoff
#' @return numerically calculated E(t) over the time vector, a matrix with rows representing time and columns representing fitness types.
#' @export
integrate_E <- function(f_initial, T_vector, argument, d_t, non_negativity_cutoff) 
{
  x_1 <- length(T_vector)     # x_1 is the dimension of the time vector
  x_2 <- length(f_initial)    # x_2 is the dimension of the type (fitness state)
  sol <- matrix(0,x_1,x_2)    # x_1 by x_2 matrix with all zeros, row is the time, column is the type (fitness state)
  sol[1,] <- f_initial        # at time zero, equals to (1-sampling probability)
  f <- f_initial              # initial E at time zero, a vector of non-sampled probabilities with the length of the types
  t=T_vector[1]               # start from time zero (current or leaves)
  
  for (ti in 2:x_1) {         # for each time point, update the vector of E
    tnext=T_vector[ti]        # next time point in the time vector
    while (t < tnext){        # update for each small time increment
      h <- min(d_t, tnext-t)  # step size of time
      k1 <- derivative_E(f, argument)
      k2 <- derivative_E(f+0.5*h*k1, argument)
      k3 <- derivative_E(f+0.5*h*k2, argument)
      k4 <- derivative_E(f+h*k3, argument)
      t <- t+h
      f <- f+h/6*(k1+2*k2+2*k3+k4)
      f[f<non_negativity_cutoff] <- non_negativity_cutoff
      #message(paste(f, collapse=" "))
    }
    sol[ti,] <- f
  }
  return (sol) 
}



#' calculate the derivative of Extinction density (a linear mutation model, acquire mutations in order)
#'
#' @param phi value of the function E 
#' @param argument a list of birth, death and mutation data
#' @return vector of calculated derivatives of E(t) 
#' @export
derivative_E <- function(phi, argument)
{
  b <- argument[[1]]
  d <- argument[[2]]
  mu <- argument[[3]]
  x_1 <- length(phi)
  dp <- numeric(x_1) # x_1 dimension vector of zeros

  #vectorized calculation for Ei (i: index of types) 
  dp[1:(x_1-1)] <- -(b[1:(x_1-1)]+d[1:(x_1-1)]+mu[1:(x_1-1)])*phi[1:(x_1-1)] + d[1:(x_1-1)] + b[1:(x_1-1)]*phi[1:(x_1-1)]*phi[1:(x_1-1)] + mu[1:(x_1-1)]*phi[1:(x_1-1)]*phi[2:(x_1)]
  dp[x_1] <- -(b[x_1]+d[x_1])*phi[x_1] + d[x_1] + (b[x_1])*phi[x_1]*phi[x_1]    # no mutation for the last state
  return (dp)
}



#' call the function integrate_E to integrate E(t)
#'
#' @param rho sampling probability
#' @param T_vector time discretization
#' @param argument a list of birth, death and mutation data
#' @param d_t step size for integration
#' @param non_negativity_cutoff
#' @return numerically calculated solution of E(t)
#' @export
integrate_phi_E <- function(rho, T_vector, argument, d_t, non_negativity_cutoff) {
  
  b <- argument[[1]]
  x_1 <- length(b)
  phi0 <- rep(1,x_1)*(1-rho)       # a vector of non-sampled probability with the length of the types
  
  sol <- integrate_E(phi0, T_vector, argument, d_t, non_negativity_cutoff) 
  
  return (sol)  
}




#' integrate_D is a function to solve differential equations using Runge-Kutta methods (https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
#' D(t) can be viewed as the density that a node at state i leads to the sampled (sub-)tree
#' similar as the "branch propagator" defined in Neher's paper (given an ancestor with fitness y at time t,
#' the prob that the child is with fitness x and time t') 
#'
#' @param f_initial vector of initial values (for each fitness type)
#' @param T_vector time discretization
#' @param argument a list of birth, death and mutation data
#' @param E_list precalculated E(t), across a list of different types
#' @param ini_rank only keeps the possible transition between state
#' @param t_1 time of the node closer to the leaves ("down node"), this is to offset E(t) which is dependent on time
#' @param d_t step size for integration
#' @param non_negativity_cutoff
#' @return numerically calculated D(t) over the time vector, a matrix with rows representing time and columns representing fitness types.
#' @export
integrate_D <- function(f_initial, T_vector, argument, E_list, ini_rank, t_1, d_t, non_negativity_cutoff) 
{
  b <- argument[[1]]                 # vector of birth rates of types
  fitness_count <- length(b)         # number of types
  x_1 <- length(T_vector)            # x_1 is the dimension of the time vector
  x_2 <- length(f_initial)           # x_2 is the dimension of the type (fitness state)  
  sol <- matrix(0,x_1,x_2)           # x_1 by x_2 matrix with all zeros, row is the time, column is the type (fitness state)
  sol[1,] <- f_initial               # at time zero, equals to f_initial
  f <- f_initial                     # initial D at time zero
  t=T_vector[1]                      # start from time zero

  for (ti in 2:x_1) {
    tnext=T_vector[ti]
    
    
    ## this step can be refined
    
    while (t < tnext){
      h <- min(d_t, tnext-t)
      
      eta_1=replicate(fitness_count, 0)  
      eta_2=replicate(fitness_count, 0)
      eta_3=replicate(fitness_count, 0)
      for (j in 1:fitness_count) {
        eta_1[j]=E_list[[j]](t+t_1)
        eta_2[j]=E_list[[j]](t+h/2+t_1)
        eta_3[j]=E_list[[j]](t+h+t_1)
      }
      
      k1 <- derivative_D(f, eta_1, argument, ini_rank)
      k2 <- derivative_D(f+0.5*h*k1, eta_2, argument, ini_rank)
      k3 <- derivative_D(f+0.5*h*k2, eta_2, argument, ini_rank)
      k4 <- derivative_D(f+h*k3, eta_3, argument, ini_rank)
      t <- t+h
      f <- f+h/6*(k1+2*k2+2*k3+k4)
      f[f<non_negativity_cutoff] <- non_negativity_cutoff
    }
    sol[ti,] <- f
  }
  return (sol)  
}




#' calculate the derivative of D
#'
#' @param phi value of the function D
#' @param eta pre-calculated Ei
#' @param argument a list of birth, death and mutation data
#' @param ini_rank only keeps the possible transition between states (the density from higher rank to lower rank is zero)
#' @return vector of calculated derivatives of D
#' @export
derivative_D <- function(phi, eta, argument, ini_rank) 
{
  b <- argument[[1]]         # birth rates
  d <- argument[[2]]         # death rates
  mu <- argument[[3]]        # mutation rates
  x_1 <- length(phi)         # number of the types
  dp <- numeric(x_1)         # x_1 dimension vector of zeros (number of types)
  dp[1:(x_1-1)] <- ini_rank[1:(x_1-1)]*(-((b[1:(x_1-1)]+d[1:(x_1-1)]+mu[1:(x_1-1)])*phi[1:(x_1-1)]) + 2*b[1:(x_1-1)]*phi[1:(x_1-1)]*eta[1:(x_1-1)] + mu[1:(x_1-1)]*eta[2:x_1]*phi[1:(x_1-1)] + mu[1:(x_1-1)]*eta[1:(x_1-1)]*phi[2:x_1])
  dp[x_1] <- ini_rank[x_1]*(-(b[x_1]+d[x_1])*phi[x_1] + 2*(b[x_1])*phi[x_1]*eta[x_1]) # no mutation for the last state
  return (dp)  
}




#' call the function integrate_D to integrate D rootward staring at time t_1 
#'
#' @param rho sampling probability
#' @param T_vector time discretization
#' @param argument a list of birth, death and mutation data
#' @param E_list precalculated E_list (see the next step)
#' @param ini_rank only keeps the possible transition between state
#' @param f_initial vector of initial values (for each type of cells)
#' @param t_1 time of the down node (starting node)
#' @param d_t step size for integration
#' @param non_negativity_cutoff
#' @return numerically calculated solution of D(t)
#' @export
integrate_phi_D <- function(rho, T_vector, argument, E_list, ini_rank, f_initial, t_1, d_t, non_negativity_cutoff){
  
  b <- argument[[1]]      # vector of birth rates of types
  x_1 <- length(b)        # number of types
  
  sol <- integrate_D(f_initial, T_vector, argument, E_list, ini_rank, t_1, d_t, non_negativity_cutoff) 
  
  return (sol)
}



#' prepare the functions of D(t) for each type transition starting from a given time t_1 
#'
#' @param rho sampling probability
#' @param T_vector time discretization
#' @param d_t step size for integration
#' @param argument a list of birth, death and mutation data
#' @param E_list precalculated E_list
#' @param t_1 time of the down node (starting node)
#' @param non_negativity_cutoff
#' @return Interpolated functions of D(t) encapsulated in a two dimensional list. the first dimension being the down node (child node) and the second being the upper node (ancestral node).
#' @export
get_D_list <- function(rho, T_vector, d_t, argument, E_list, t_1, non_negativity_cutoff) {
  b <- argument[[1]]                      # vector of birth rates of types
  fitness_count <- length(b)              # number of types
  D_list <- list()                        # list of functions D_i(t)
  
  for (i in 1:fitness_count) {            # for each fitness type (CHILD)
      
    phi0 <- rep(0, times = fitness_count) # start all other types with 0
    if (t_1==0){                          # at a leaf, set the initial value to be the sampling probability
      phi0[i] <- rho                      
    } else {                              # at an internal node, set the initial value to be 1
      phi0[i] <- 1
    }
    

    # the ini_rank is needed because a cell can only mutate into the higher fitness state, NO return to the lower fitness state
    # as such, type 2 CANNOT mutate to type 1
    temp1 <- rep(1, times = i)                      ## initialize ini_rank, 1s for up to the current type
    temp2 <- rep(0, times = (fitness_count-i))      ## initialize ini_rank, 0s for the higher states
    ini_rank <- c(temp1, temp2)                     ## initialize ini_rank, combine

    D_sol <- integrate_phi_D(rho, T_vector, argument, E_list, ini_rank, phi0, t_1, d_t, non_negativity_cutoff)
    temp_list <- list()
    for (j in 1:fitness_count) {                                    # for each fitness type (ANCESTRAL)
      #if (sum(is.na(D_sol[,j]))>100){
        #browser()
      #}
      D_approx <- approxfun(T_vector, D_sol[,j], method="linear")   # approximate the functions of D_i(t) by interpolation
      temp_list[[j]] <- D_approx
    }
    D_list[[i]] <- temp_list   # i represents the type of the node closer to the leaves
  }
  return (D_list)  
}




#' calculate the propagator for a branch, a propagator is the probability density from an ancestral node to a child
#'
#' @param rho sampling probability
#' @param argument a list of birth, death and mutation data
#' @param t time between up node and down node (length of the branch)
#' @param t_1 starting time from the child
#' @param E_list pre-calculated E_list (see the next step)
#' @param T_vector time discretization
#' @param d_t step size for integration
#' @param non_negativity_cutoff
#' @return a matrix of integrated D values of a branch (ancestral type to the child type)
#' @export
integrate_prop <- function(rho, argument, t, t_1, E_list, T_vector, d_t, non_negativity_cutoff){
  b <- argument[[1]]
  fitness_count <- length(b)
  
  D_list <- get_D_list(rho, T_vector, d_t, argument, E_list, t_1, non_negativity_cutoff)
  
  sol <- replicate(fitness_count, numeric(fitness_count))
  
  for (i in 1:fitness_count){
    for (j in 1:fitness_count){
      sol[i,j] <- D_list[[j]][[i]](t)
    }
  }
  return (sol)  
}


