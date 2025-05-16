clear;clc

%This script imports and processes the simulation data from the P1C vs. 
%P1P2 scan

sim_name = 'P1C_vs_P1P2_linear_sweep';
num_jobs = 120;
run_per_job = 10;
E_cutoff = 1e-6;
enzyme_only = 1;
[parameter_table,sr] = import_batch_run(sim_name,num_jobs,run_per_job,E_cutoff);

sr.burden_mat = zeros(size(sr.raw_c_cell));

for l = 1:size(sr.raw_c_cell,1)
    for m = 1:size(sr.raw_c_cell,2)
        rho = sr.raw_rho_cell{l,m};
        c = sr.raw_c_cell{l,m};
        burden = 0;
        for n = 1:length(rho)
            burden = burden + rho(n)*sum(c(n,:));
        end
        sr.burden_mat(l,m) = burden/sum(rho);    
    end
end
        
save(['../../../figure_generation/fig_data/',sim_name,'_simulation_results.mat'])