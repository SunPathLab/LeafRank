

    tree_name = "DM-tree.mat"
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
    r_ub = 2.5; 
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
    
    
    
    
    
    
    %% Optimization
    
    options_g = optimoptions('fmincon','Algorithm','interior-point', 'MaxIterations',1e4, 'MaxFunctionEvaluations',1e5,'Display','off','SpecifyObjectiveGradient',true);
    
    fun_g_WGD = @(x) poisson_WGD_grad(x,edges,li,ni,edge_wgd,edge_length);
    f_hist = [];
    x_hist = [];
    for ii = 1:N
        try
            [x_g,f_g] = fmincon(fun_g_WGD, x_init(:,ii), A, b, [], [], LOW, UP, [], options_g);
            f_hist = [f_hist,f_g];
            x_hist = [x_hist,x_g];
        catch
        end
    end
    
    [mf,idx] = min(f_hist);
    
    
    writematrix(x_hist(:,idx),"ultra-tree.csv")
    
    save("ultra-tree.mat")



