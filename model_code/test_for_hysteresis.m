function [ET_store,simu_store] = ...
    test_for_hysteresis(environments,table_file)

%This function performs a test for hysteresis across specified
%environments. Rather than starting from an abiotic initial state, it
%starts from the end state of other simulations. 

global par

%Load parameter table
par_table = readtable(table_file,'ReadRowNames',false,...
    'ReadVariableNames',true);

%Replace parameters
non_id_vars = par_table.Properties.VariableNames(1:(end-1));
matching_row = par_table.id == environments(1);
for i = 1:length(non_id_vars)
    par.(non_id_vars{i}) = par_table{matching_row,non_id_vars{i}};
end

%Set up storage vectors
ET_store = cell(size(environments));
simu_store = ET_store;

%Simulate first environment
[ET_store{1},~, ~,~,~,simu_store{1}] = chemostat_invasion_ode15s();

%Loop through and simulate remaining environments 
par.existing_environment = true;
for i = 2:length(environments)
    
    %Set up parameters for this simulation
    non_id_vars = par_table.Properties.VariableNames(1:(end-1));
    matching_row = par_table.id == environments(i);
    for j = 1:length(non_id_vars)
        par.(non_id_vars{j}) = par_table{matching_row,non_id_vars{j}};
    end

    %Simulate steady state with existing species but new conditions
    [~,~,~,~,~,simu] = simulate_crossfeeding_ode15s(simu_store{i-1});
    par.existing_rho = simu.rho;
    par.existing_c = simu.c;
    
    %Perform evolutionary invasions to find new ESS
    [ET_store{i},~,~,~,~,simu_store{i}] = chemostat_invasion_ode15s();

end

end

