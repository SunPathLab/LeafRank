function mutation_data=find_mutation(mutation_file, ancestry_data)

% generate a mutation file matrix for sampled cells

[row_1, col_1]=size(ancestry_data);
[row_2, ~]=size(mutation_file);

mutation_data=[];


index_vec=zeros(row_1,1);

for i=1:row_2  % iterate through mutation file
    a=mutation_file(row_2+1-i, 1);
    for j=1:row_1
        for k=2:col_1
            if ancestry_data(j,k)==a
                index_vec(j)=index_vec(j)+1;
                mutation_data(j,index_vec(j))=mutation_file(row_2+1-i, 2);
            end
        end
    end
end
        