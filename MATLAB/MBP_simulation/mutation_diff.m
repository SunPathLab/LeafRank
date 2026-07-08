function distance_matrix=mutation_diff(mutation_data)

[row_1, col_1]=size(mutation_data);
temp=zeros(row_1, col_1+1);
temp(1:row_1, 1:col_1)=mutation_data;


distance_matrix=zeros(row_1, row_1);


for i=1:row_1
    for j=1:row_1
        diff_1=setdiff(temp(i,:),temp(j,:));
        diff_2=setdiff(temp(j,:),temp(i,:));
        distance_matrix(i,j)=size(diff_1,2)+size(diff_2,2);
    end
end
        