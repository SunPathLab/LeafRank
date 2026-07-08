function [ancestry_data, node_1, node_2, node_3, edge_length_1, edge_length_2, node_index, current_high, pop_index]=find_merge(ancestry_data, node_index, current_high, pop_index, time_history, T)

% edited 2025/06/03

%% get the most recent merge information from the ancestry data

% node_index: index of nodes for each cell in the tree




[~, col_1]=size(ancestry_data);

is_done=0;

while is_done<1

    maximum = max(ancestry_data(:,1));
    max_loc=find(ancestry_data(:,1)==maximum);
    [row_temp, ~]=size(max_loc);

    if row_temp==1 %% if no merge happens at the cell with the maximum index
        temp_data=ancestry_data(max_loc,:);
        temp_data=temp_data(2:col_1);
        temp_data(col_1)=0;
        ancestry_data(max_loc,:)=temp_data;
    else
        % a merge happens
        is_done=1;

        node_1=node_index(max_loc(1));
        node_2=node_index(max_loc(2));

        % update node_index
        node_index(max_loc(2))=[];
        node_index(max_loc(1))=[];
        current_high=current_high+1;
        node_index=[node_index, current_high];

        node_3=current_high;

        pop_index_1=pop_index(max_loc(1));
        pop_index_2=pop_index(max_loc(2));
        pop_index_3=maximum;

        if time_history(pop_index_1, 2)==0 % the cell is still alive
            edge_length_1=T-time_history(pop_index_3, 2);
        else
            edge_length_1=time_history(pop_index_1, 2)-time_history(pop_index_3, 2);
        end

        if time_history(pop_index_2, 2)==0
            edge_length_2=T-time_history(pop_index_3, 2);
        else
            edge_length_2=time_history(pop_index_2, 2)-time_history(pop_index_3, 2);
        end

        temp_data=ancestry_data(max_loc(2),:);
        ancestry_data(max_loc(2),:)=[];
        ancestry_data(max_loc(1),:)=[];
        ancestry_data=[ancestry_data; temp_data];

        pop_index(max_loc(2),:)=[];
        pop_index(max_loc(1),:)=[];
        pop_index=[pop_index; maximum];
    end
end

