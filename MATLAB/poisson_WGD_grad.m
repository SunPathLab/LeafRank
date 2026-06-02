function [obj,grad] = poisson_WGD_grad(param, edges, li, ni, wgd_idx, branch_length)
    obj     = poisson_WGD(param,edges, li, wgd_idx, branch_length);

    rate    = param(1:2);
    time = [zeros(li,1);1;param(3:end)];
    r_time  = time(edges(:,1)) - time(edges(:,2));
    grad    = zeros(length(param),1);
    r1_idx  = find(wgd_idx == 0);
    r2_idx  = find(wgd_idx == 1);
    grad(1) = sum(branch_length(r1_idx)./rate(1) - r_time(r1_idx));
    grad(2) = sum(branch_length(r2_idx)./rate(2) - r_time(r2_idx));
    
    for i = 1:ni-1
        node_id = li + 1 + i;
        child   = find(edges(:,1) == node_id);
        c1 = child(1);
        c2 = child(2);
        p  = find(edges(:,2) == node_id);
        if length(p) > 1
            keyboard
        end
        grad(2 + i) = rate(wgd_idx(p)+1) - rate(wgd_idx(c1)+1) - rate(wgd_idx(c2)+1) + branch_length(c1)/r_time(c1) + branch_length(c2)/r_time(c2) - branch_length(p)/r_time(p);
    end
    grad = -grad;
end