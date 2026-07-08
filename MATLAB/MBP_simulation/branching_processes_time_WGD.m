function [time, type]=branching_processes_time_WGD(birth_rate, death_rate, state, WGD_state, passenger_rate, WGD_driver)

    % edited 2025/11/05
    
    % birth_rate: a matrix of birth rates
    % death_rate: a vector of death rates
    % passenger_rate: rate of passenger mutation occurred during the lifetime of a cell
    % state: state of the cell (or type of the cell)
    % type_count: the number of cell states
    % WGD_driver: coefficient for WGD event driver effect
    


    type=zeros(1,4);
    if WGD_state >= 2
        birth_rate = 3*birth_rate - 2*diag(diag(birth_rate));
    end
    % the state of target cell, the state of child, event type, WGD state
    % event type 0: death event
    % event type 1: birth event
    % event type 2: passenger mutation event during the lifetime
    
    total_rate=0;
    for i=1:size(birth_rate,2)
        total_rate=total_rate+birth_rate(state, i);
    end
    total_rate=total_rate+death_rate(state);
    for i=1:size(passenger_rate,2)
        try
            total_rate=total_rate+passenger_rate(WGD_state,i);
        catch
            keyboard
        end
    end
    % obtain the total rate
    time=exprnd(1/total_rate);
    
    u=rand;
    
    check=0;
    temp_1=0;
    
    
    for i=1:size(birth_rate,2)
        temp_1=temp_1+birth_rate(state, i)/total_rate;
        if u<temp_1
            check=1;
            type=[state, i, 1, WGD_state];  % birth event, gives birth to a type i cell
            break
        end
    end
    
    if check==0
        temp_1=temp_1+death_rate(state)/total_rate;
        if u<temp_1
            check=1;
            type=[state, state, 0, WGD_state];  % death event
        end
    end
    
    for i=1:size(passenger_rate,2)
        if check==0
            
            temp_1 = temp_1 + passenger_rate(WGD_state,i)/total_rate;
            if u < temp_1
                check=1;
                if WGD_state == 1 && i == 2
                    type = [state+WGD_driver, state+WGD_driver, 2, i]; % WGD event
                else
                    type=[state, state, 2, i];  % passenger mutation event
                end
            end
        end
    end
    


end

    
    




