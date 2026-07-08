function LL = poisson_WGD(param, edges, li, wgd_idx, branch_length)
    rate = param(1:2);
    time = [zeros(li,1);1;param(3:end)];
    r_time = time(edges(:,1)) - time(edges(:,2));
    % if len(time) ~= len(wgd_idx)
    %     keyboard
    % end
    Xi = r_time.*rate(wgd_idx + 1);
    LL = -sum(branch_length.*log(Xi) - Xi );
end