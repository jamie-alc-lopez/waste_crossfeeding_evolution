function [transition_points,transition_map,current_bounds] = ...
    find_linear_phase_boundaries(S,n_points,n_bisections,alpha1,alpha2,...
    beta,D_range,Vc,Keq1,Keq2,gm)

%This script computes phase boundaries based on the analytical approximation
%of the linear thermodynamics model 

%Define bisection function
bfun = @(x) compute_ESS_consortium(S,alpha1,alpha2,beta,x,Vc,Keq1,Keq2,gm);

%Get initial supply vector
D_vec = linspace(D_range(1),D_range(2),n_points);

%Get initial ESS consortia
for i = 1:n_points
    [ESS_ind(i),ESS_type{i},optimal_ET{i},optimal_c1e{i},ss_vals{i}] = bfun(D_vec(i));
end

%Identify preliminary transition points
ind_diff = diff(ESS_ind);
transition_ind = find(ind_diff);
n_transitions_init = length(transition_ind);
transition_map = nan(n_transitions_init,2);
current_bounds = nan(n_transitions_init,2);

%Make transition map single index if no transition
if n_transitions_init == 0
    transition_map = [ESS_ind(1),ESS_ind(1)];
    current_bounds = [0,0];
end

%Loop through transitions
loop_var = n_transitions_init > 0;
trn_id = 1;
while loop_var
    
    %If transition was id'd initially, get info
    if trn_id <= n_transitions_init
        ind_set = [transition_ind(trn_id),transition_ind(trn_id)+1];
        transition_map(trn_id,:) = [ESS_ind(ind_set(1)),ESS_ind(ind_set(2))];
        current_bounds(trn_id,:) = [D_vec(ind_set(1)),D_vec(ind_set(2))];
    end
    
    %Get updated
    [current_bounds,transition_map] = ...
        bisect_transition(n_bisections,bfun,trn_id,transition_map,current_bounds);
    
    %Update number of transitions
    n_transitions = size(transition_map,1);
    
    %Check if loop exit criteria met, else advance counter
    if trn_id >= n_transitions
        loop_var = 0;
    else
        trn_id = trn_id + 1;
    end
    
end

%Get middle points (true transition points)
transition_points = mean(current_bounds,2);

end

