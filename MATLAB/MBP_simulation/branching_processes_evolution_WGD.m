function [population, z_0, mutation_num_passenger, cell_tot, parent, mutation_file, T, time_history, WGD_history]=branching_processes_evolution_WGD(population, z_0, birth_rate, death_rate, mutation_num_passenger, cell_tot, parent, mutation_file, passenger_rate, time_history, WGD_history, WGD_driver)

    % edited 2025/06/03
    
    %% population: stores information of all living cells
    %% z_0: total number of living cells
    
    %% population is ranked in time of the next event
    p_temp=population.peek();
    
    t=p_temp(1,1); % time of the next event
    
    type=p_temp(1,2:5);
    % the state of target cell, the state of child, event type
    % event type 0: death event
    % event type 1: birth event
    % event type 2: passenger mutation event during the lifetime
    
    parent_id=p_temp(1,6); % cell index of the cell responsible for the next event
    
    T=t; % update time
    
    
    time_history(parent_id, 2)=T;
    % time_history stores the birth time and death time of each cell
    % update the "death time"
    
    population.remove(); % remove the cell at top
    
    
    
    if type(3)==0 % death happens
        z_0=z_0-1;
     
    
    elseif type(3)==1 % birth happens
        state=type(1);
        WGD_state=type(4);
        cell_tot=cell_tot+1;
        [time, type_1]=branching_processes_time_WGD(birth_rate, death_rate, state,WGD_state, passenger_rate, WGD_driver);
        temp_1=[time+T, type_1, cell_tot];
        
        time_history(cell_tot, 1)=T;
        
        population.insert(temp_1);
        parent(end+1,1)=parent_id;
    
    
        state=type(2);
        cell_tot=cell_tot+1;
        [time, type_2]=branching_processes_time_WGD(birth_rate, death_rate, state,WGD_state, passenger_rate, WGD_driver);
        temp_2=[time+T, type_2, cell_tot];
        
        time_history(cell_tot, 1)=T;
    
        population.insert(temp_2);
        parent(end+1,1)=parent_id;
    
        z_0=z_0+1;
       
    else % passenger mutation during lifetime
        state=type(1);
        WGD_state=type(4);
        cell_tot=cell_tot+1;
        [time, type_1]=branching_processes_time_WGD(birth_rate, death_rate, state,WGD_state, passenger_rate, WGD_driver);
        temp_1=[time+T, type_1, cell_tot];
        
        time_history(cell_tot, 1)=T;
    
        population.insert(temp_1);
        parent(end+1,1)=parent_id;
    
        mutation_num_passenger=mutation_num_passenger+1;
        mutation_file(end+1,:)=[cell_tot, mutation_num_passenger];
        WGD_history(end+1,:)=[cell_tot, type_1(4)];
    end

end

               



 


            
            
            
            
    