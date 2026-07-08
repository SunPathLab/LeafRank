function MBP_simulation(output_dir, options)
    %% Extract the simulation setting 
    dir_name = output_dir;
    if ~isfolder(dir_name)
        mkdir(dir_name)
    end

    
    if isfield(options, "tree_sample_size")
        tree_sample_size = options.tree_sample_size;
    else
        tree_sample_size=250;
    end

    if isfield(options, "passenger_rate")
        passenger_rate = options.passenger_rate;
    else
        passenger_rate = 0.18;  
    end

    if isfield(options, "WGD_rate")
        WGD_rate = options.WGD_rate;
    else
        WGD_rate = 0;  
    end

    if isfield(options, "num_types")
        num_types = options.num_types;
    else
        num_types = 8;
    end
    
    if isfield(options, "driver_rate")
        driver_rate = options.driver_rate;
    else
        driver_rate = 0.0001;
    end

    if isfield(options, "birth_rates")
        birth_rates = options.birth_rates;
        if length(birth_rates) ~= number_types
            error("Number of birth rates should be the same as the number of types")
        end
    else
        birth_rates = linspace(0.2,0.8,num_types);
    end

    if isfield(options, "death_rate")
        death_rate = options.death_rate;
    else
        death_rate = 0.18;
    end

    if isfield(options, "total_cells")
        n = options.total_cells;
    else
        n = 1000000;
    end

    

    %% simulate the tumor
    
    WGD_driver = 0;  % WGD driver fitness effect index. This variable is designed for evaluating simulation under potential fitness effect of WGD. Therefore it should be set as 0 for current publication.
    driver_rates = driver_rate*ones(num_types-1,1);
    birth_mat = diag(birth_rates) + diag(driver_rates,1);
    death_rates = death_rate*ones(num_types,1);
    aberration_mat = [passenger_rate, WGD_rate;0,2*passenger_rate]; % Here WGD is assumed to double the passenger mutation rates.

    [population, z_0, T, ~, ~, parent, mutation_file, time_history, WGD_history]=branching_processes_main_WGD(birth_mat, death_rates, n, aberration_mat, WGD_driver);
    temp_tree=randsample(n,tree_sample_size);
    
    %%
    
    
    p_info=zeros(n,6);
    for j=1:n
        p_info(j,:)=population.peek();
        population.remove();
    end
    
    
    sample_tree=p_info(temp_tree,6);
    ancestry_data=find_ancestry(parent, sample_tree);
    mutation_data=find_mutation(mutation_file, ancestry_data);
    [row_1,col_1]=size(mutation_data);
    mutation_data(row_1+1,1:col_1+1)=zeros(1,col_1+1);
    distance_matrix=mutation_diff(mutation_data);
    str_1 = strcat(dir_name,'/DM.txt');
    writematrix(distance_matrix, str_1);
    
    
    str_2 = strcat(dir_name,'/FT.txt');
    fitness=p_info(temp_tree,2);
    writematrix(fitness, str_2);
    
    edge_matrix=[];
    length_matrix=[];
    [row_2,~]=size(ancestry_data);
    node_index=[1:1:row_2]; %% node index for leaves
    current_high=row_2; %% store the highest node index so far
    pop_index=sample_tree;
    for k=1:(current_high-1)
        [ancestry_data, node_1, node_2, node_3, edge_length_1, edge_length_2, node_index, current_high, pop_index]=find_merge(ancestry_data, node_index, current_high, pop_index, time_history, T);
        edge_matrix=[edge_matrix; node_3, node_1];
        edge_matrix=[edge_matrix; node_3, node_2];
        length_matrix=[length_matrix; edge_length_1];
        length_matrix=[length_matrix; edge_length_2];
    end
    edge_matrix=[edge_matrix; current_high+1, current_high];
    length_matrix=[length_matrix; time_history(ancestry_data(1, 1),2)];
    str_3 = strcat(dir_name,'/EDGE.txt');
    writematrix(edge_matrix, str_3);
    str_4 = strcat(dir_name,'/LEN.txt');
    writematrix(length_matrix, str_4);
    str_5 = strcat(dir_name,'/WGD.txt');
    writematrix(p_info(temp_tree,5),str_5);
    
    
    save(strcat(dir_name,'/data.mat'))
end