# MATLAB Helper Scripts for R package LeafRank

This directory contains two major MATLAB-based helper modules that support the analysis of LeafRank:
1. Non-spatial simulation model `MBP_simulation`:
    This directory contains all functions required to simulate the non-spartial synthesis tumors based on a multi-type branching process model. The simulation is performed stochastically using the Gillespie algorithm, which randomly generates the elapsed time between cellular events. These events are stored in a priority queue to preserve their chronological order. A detailed description of the modeled events can be found in the manuscript

2. Modified chronos implementation `modified_chronos`:
    This directory contains all functions required to convert a distance-matrix based tree into an ultrametric tree under the state-dependent molecular clock model proposed in the manuscript. It requires inputs contains
    - node_idx: index for internal nodes
    - leaf_idx: index for leaf nodes
    - edge_length: distance of each edge
    - edge_wgd: indicator for edge WGD status
    - edges: edges information presented by ancestor and descendent nodes index


We provide two MATLAB script, `MBP_simulation_example.m` and `modified_chronos_example.m`, to demonstrate the usage of these two modules. In each script, users can specify the following options and then run the entire script:

1. `MBP_simulation_example.m`:
    - tree_sample_size: Scalar specifying the number of cells sampled to construct the phylogenetic tree.
    - passenger_rate: Scalar specifying the rate of passenger aberrations.
    - WGD_rate: Scalar specifying the rate of whole-genome doubling. This can be set to 0 to simulate a non-WGD scenario.
    - num_types: Integer scalar specifying the number of potential fitness types, each characterized by a distinct cellular birth rate.
    - driver_rate: Scalar specifying the rate of driver aberrations. This is assumed to be identical across all fitness types.
    - birth_rates: 1 x num_types vector specifying the birth rates of the distinct fitness types.
    - death_rate: Scalar specifying the cellular death rate. This is assumed to be identical acorss all fitness types.
    - total_cells: Scalar specifying the number of simulated live cells. The simulation stops once the number of live cells reaches this value.

2. `modified_chronos_example.m`:
    - rate_bounds: 1 x 2 vector specifying the lower and upper bounds for aberration rates.
    - WGD_bounds: 1 x 2 vector specifying the lower and upper bounds for the multiplier applied to aberration rates after whole-genome doubling.
    - init_num: Scalar specifying the number of initial points used for optimization. Due to non-convexity, a multi-start approach is implemented to identify the global optimum.
    - opt_algorithm: String specifying the optimization solver used in the `fmincon` nonlinear constrained optimization function. This requires the MATLAB `Optimization Toolbox`.
    - opt_max_iter: Scalar specifying the maximum number of iterations in the optimization function.
    - is_par: Boolean indicator for parallelization. Parallelization is recommended due to the complexity of the optimization procedure. This requires the MATLAB `Parallel Computing Toolbox`.


These modules were developed and tested using MATLAB version 24.2.0.2712019 (R2024b) with the following add-on toolboxes: `Optimization Toolbox` and `Parallel Computing Toolbox`