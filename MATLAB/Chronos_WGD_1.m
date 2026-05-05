%case_num  = '8';
%batch_num = '42';
%dir_name  = strcat("Results/SIM-spatial-",batch_num,"-250/run-1");
OV_num    = "075";
dir_name  = strcat("Results/Subclonal/OV-",OV_num);
%case_num  = '42';
parpool('local',50)
warning('off','all')

%for i = 25:50

    %batch_num = num2str(i);
    %dir_name  = strcat('Results/SIM-spatial-',case_num,'/run-',batch_num);
    tree_name = strcat(dir_name,"/tree_alt_dec.mat");
    load(tree_name)
    
    
    %% Setting the constraints
    li = length(leaf_idx);
    ni = length(node_idx);
    ages = zeros(li+ni,1);
    ages(li +1) = 1;
    
    rate_lb = 1e-6;
    rate_ub = 1e5;
    time_lb = 1e-8;
    time_ub = 1;
    
    LOW = [rate_lb;rate_lb;time_lb*ones(ni-1,1)];
    UP  = [rate_ub;rate_ub;time_ub*ones(ni-1,1)];
    
    
    A = zeros(2,length(LOW));
    b = zeros(2,1);
    r_lb = 1.03;
    r_ub = 5;
    A(1:2,1:2) = [r_lb,-1;-r_ub,1];
    
    for e = 1:size(edges,1)
        temp = zeros(1,ni+1);
        if edges(e,1) <=li+1 || edges(e,2) <= li+1
            continue
        end
        temp(edges(e,2)-li+1) = 1;
        temp(edges(e,1)-li+1) = -1;
        A = [A;temp];
        b = [b;0];
    end
    
    
    
    
    
    
    %% Randomized the initial times
    known_ages = linspace(1,li,li);
    
    N  = 100;
    ages = zeros(li+ni,N);
    
    t0 = zeros(ni,N);
    
    while length(known_ages) < li+ni
        for i = 1:length(t0)
            if t0(i,:) ~= 0
                continue
            end
            c = find(edges(:,1) == li+i);
            c1 = edges(c(1),2);
            c2 = edges(c(2),2);
            if ismember(c1,known_ages) && ismember(c2,known_ages) 
    
                ages(li+i,:) = max(ages(c1),ages(c2)) + rand(1,N);
                t0(i,:) = ages(li+i);
                known_ages = [known_ages,li+i];      
            end
    
        end
    end
    
    t0 = t0./max(t0);
    
    x_init = [];
    for i = 1:N
        x1 = LOW(1) + rand*(UP(1)-LOW(1));
        x2 = (r_lb + rand*(r_ub-r_lb))*x1;
        x_init = [x_init,[x1;x2;t0(2:end,i)]];
    end
    
    
    
    
    
    
    %%
    
    %options = optimoptions('fmincon','Algorithm','sqp', 'MaxIterations',1e5, 'MaxFunctionEvaluations',1e5);
    options_g = optimoptions('fmincon','Algorithm','interior-point', 'MaxIterations',1e4, 'MaxFunctionEvaluations',1e5,'Display','off','SpecifyObjectiveGradient',true);
    
    fun_WGD = @(x) poisson_WGD(x,edges, edge_wgd, edge_length);
    fun_g_WGD = @(x) poisson_WGD_grad(x,edges,li,ni,edge_wgd,edge_length);
    %[x,f] = fmincon(fun_WGD, x0,A,b,[],[],LOW,UP,[],options);
    f_hist = [];
    x_hist = [];
    parfor ii = 1:N
        try
            [x_g,f_g] = fmincon(fun_g_WGD, x_init(:,ii), A, b, [], [], LOW, UP, [], options_g);
            f_hist = [f_hist,f_g];
            x_hist = [x_hist,x_g];
        catch
        end
    end
    
    [mf,idx] = min(f_hist);
    
    
    writematrix(x_hist(:,idx),strcat(dir_name,"/ultra-tree_alt_dec_100.csv"))
    
    save(strcat(dir_name,"/ultra-tree_alt_dec_100.mat"))
    % sol_hist = [];
    % f_hist   = [];
    % for i = 1:10
    %     [x,f] = fmincon(fun_WGD, x_init(:,i),A,b,[],[],LOW,UP,[],options);
    %     sol_hist = [sol_hist,x];
    %     f_hist = [f_hist,f];
    % end
    
    
    %%
    % parents = edges(:,1);
    % children = edges(:,2);
    % 
    % root = setdiff(parents, children);   % root node(s)
    % leaves = setdiff(children, parents); % leaf nodes
    % 
    % paths = cell(numel(leaves),1);  % store paths as cell array
    % 
    % for i = 1:numel(leaves)
    %     node = leaves(i);
    %     pathEdges = [];
    % 
    %     while node ~= root
    %         % find the edge index where this node is a child
    %         idx = find(edges(:,2) == node);
    %         pathEdges = [idx; pathEdges];   % prepend edge index
    %         node = edges(idx,1);            % move to parent
    %     end
    % 
    %     paths{i} = pathEdges;   % store edge indices for this leaf
    % end

%end




