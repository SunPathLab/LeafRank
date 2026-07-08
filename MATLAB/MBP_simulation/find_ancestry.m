function ancestry_data=find_ancestry(parent, sample)

% edited 2025/06/03

%% generate an ancestry matrix for sampled cells

[row_1, ~]=size(sample);
[row_2, ~]=size(parent);
ancestry_data=sample(:,1);


index_vec=ones(row_1,1);

for i=1:row_2
    a=parent(row_2+1-i, 1); % from bottom to top
    for j=1:row_1
        if ancestry_data(j,index_vec(j))==row_2+1-i
            index_vec(j)=index_vec(j)+1;
            ancestry_data(j,index_vec(j))=a;
        end
    end
end
        
    


