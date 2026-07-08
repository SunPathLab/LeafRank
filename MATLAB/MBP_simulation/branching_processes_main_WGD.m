function [population, z_0, T, mutation_num_passenger, cell_tot, parent, mutation_file, time_history, WGD_history]=branching_processes_main_WGD(birth_rate, death_rate, n, passenger_rate, WGD_driver)

    % edited 2025/06/03
    
    % n: targeted population size
    % birth_rate: a matrix of birth rates
    % death_rate: a vector of death rates
    % passenger_rate: rate of passenger mutation
    
    % type_count: number of types of cells
    [type_count,~]=size(birth_rate);
    
    % z_0: number of living cells
    z_0=1;
    
    % cell_tot: total number of cells created in history
    cell_tot=1;
    
    % T: time of the simulation
    T=0;
    
    % mutation_num_passenger: index of passenger mutations
    mutation_num_passenger=0;
    
    % The first cell is of type 1
    [time, type]=branching_processes_time_WGD(birth_rate, death_rate, 1,1, passenger_rate, WGD_driver);
    
    population=PriorityQueue(1);
    %% initialize a priority queue to store population information
    population.insert([time, type, 1]);
    %% time of next event, type of next event (1 by 3), cell index 
    
    %% the first cell has no parent, this vector stores the index of the parent of each cell
    parent=[0]; 
    
    
    %% time_history stores the birth time and death time of each cell
    time_history=[0, 0];

    %% WGD state history stores the cell id and WGD status
    WGD_history = [1,1];
    
    
    %% mutation_file stores which cell acquires a passenger mutation
    mutation_file=[];
    
    
    % population.peek()
    % population.size()
    % population.remove()
    % iterate until the population size reaches n
    while z_0<n 
    
        [population, z_0, mutation_num_passenger, cell_tot, parent, mutation_file, T, time_history, WGD_history]=branching_processes_evolution_WGD(population, z_0, birth_rate, death_rate, mutation_num_passenger, cell_tot, parent, mutation_file, passenger_rate, time_history, WGD_history, WGD_driver);

    
          
        if z_0==0 %% restart the simulation if the population goes extinct
            cell_tot=1;
            z_0=1;
            T=0;
            mutation_num_passenger=0;
            [time, type]=branching_processes_time_WGD(birth_rate, death_rate, 1, 1, passenger_rate, WGD_driver);
            population=PriorityQueue(1);
            population.insert([time, type, 1]);
            parent=[0];
            time_history=[0, 0];
            WGD_history = [1,1];
            mutation_file=[];
        end        
    end

end