function modified_chronos(input_tree_path, output_csv_path, output_all_path, options)
    %% Normalize the Inputs 
    
    [input_dir, input_name, input_ext] = fileparts(input_tree_path);

    if input_ext == ".mat" || input_ext == ""
        input_ext = ".mat";
    else
        error("Input tree data must end with .mat");
    end

    [output_dir, output_name, output_ext] = fileparts(output_csv_path);

    if output_ext == ".csv" || output_ext == ""
        output_ext = ".csv";
    else
        error("Output csv path must end with .csv");
    end


    [output_all_dir, output_all_name, output_all_ext] = fileparts(output_all_path);

    if output_all_ext == ".mat" || output_all_ext == ""
        output_all_ext = ".mat";
    else
        error("Output all variables path must end with .mat");
    end

    if isfield(options, "rate_bounds")
        rate_lb = options.rate_bounds(1);
        rate_ub = options.rate_bounds(2);
    else
        rate_lb = 1e-6;
        rate_ub = 1e5;
    end
    

    if isfield(options, "WGD_bounds")
        WGD_lb = options.WGD_bounds(1);
        WGD_ub = options.WGD_bounds(2);
    else
        WGD_lb = 1.03;
        WGD_ub = 2.5;
    end

    if isfield(options, "init_num")
        N = options.init_num;
    else
        N = 100;
    end

    if isfield(options, "opt_algorithm")
        opt_algo = options.opt_algorithm;
    else
        opt_algo = "interior-point";
    end

    if isfield(options, "opt_max_iter")
        opt_max_iter = options.opt_max_iter;
    else
        opt_max_iter = 1e4;
    end

    if isfield(options, "is_par")
        is_par = options.is_par;
    else
        is_par = false;
    end

    
    
    options_g = optimoptions('fmincon','Algorithm',opt_algo, 'MaxIterations',opt_max_iter, 'MaxFunctionEvaluations',1e5,'Display','iter','SpecifyObjectiveGradient',true);


    tree_name = strcat(input_dir,input_name, input_ext);
    load(tree_name, 'node_idx', 'leaf_idx', 'edge_length', 'edge_wgd', 'edges')
    
    %% Setting the constraints
    li = length(leaf_idx);
    ni = length(node_idx);
    ages = zeros(li+ni,1);
    ages(li +1) = 1;
    time_lb = 1e-8;
    time_ub = 1;
    LOW = [rate_lb;rate_lb;time_lb*ones(ni-1,1)];
    UP  = [rate_ub;rate_ub;time_ub*ones(ni-1,1)];
    A = zeros(2,length(LOW));
    b = zeros(2,1);
    A(1:2,1:2) = [WGD_lb,-1;-WGD_ub,1];
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
        x2 = (WGD_lb + rand*(WGD_ub-WGD_lb))*x1;
        x_init = [x_init,[x1;x2;t0(2:end,i)]];
    end
    
    
    %% Optimization
    
    
    fun_g_WGD = @(x) poisson_WGD_grad(x,edges,li,ni,edge_wgd,edge_length);
    f_hist = [];
    x_hist = [];

    if is_par
        parfor ii = 1:N
            try
                [x_g,f_g] = fmincon(fun_g_WGD, x_init(:,ii), A, b, [], [], LOW, UP, [], options_g);
                f_hist = [f_hist,f_g];
                x_hist = [x_hist,x_g];
            catch
            end
        end
    else
        for ii = 1:N
            try
                [x_g,f_g] = fmincon(fun_g_WGD, x_init(:,ii), A, b, [], [], LOW, UP, [], options_g);
                f_hist = [f_hist,f_g];
                x_hist = [x_hist,x_g];
            catch
            end
        end
    end
    [mf,idx] = min(f_hist);
    

    
    
    writematrix(x_hist(:,idx),strcat(output_dir,output_name, output_ext))
    
    save(strcat(output_all_dir, output_all_name, output_all_ext))

end