function nonzero_sol_ind = find_nonzero_sols(sol,nonzero_vars)

%This script identifies nonzero solutions within a set of solutions
%produced by solve.

fn = fieldnames(sol);
num_sols = length(sol.(fn{1}));
num_nonzero_vars = length(nonzero_vars);
nonzero_sol_ind = ones(num_sols,1);

for i = 1:num_sols
    for j = 1:num_nonzero_vars
       is_nonzero =  ~isequaln(sol.(nonzero_vars{j})(i),sym(0));
       nonzero_sol_ind(i) = nonzero_sol_ind(i) & is_nonzero; 
    end    
end


end

