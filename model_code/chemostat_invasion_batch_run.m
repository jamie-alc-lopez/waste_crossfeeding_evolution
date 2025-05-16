function chemostat_invasion_batch_run(table_file,batch_id,results_folder,run_per_job)

% This script is a wrapper to run batches of chemostat invasion simulations
% that can be run from command line in an HPC context

global par

rng('shuffle')

%Load parameter table
par_table = readtable(table_file,'ReadRowNames',false,...
    'ReadVariableNames',true);

%Initialize storage variables
ET_cell = {};
growth_cell = {};
simu_cell = {};
num_auto_peak_cell = {};

%Compute start and end simulation indices
start_end(1) = (batch_id-1)*run_per_job + 1; 
start_end(2) = start_end(1) + run_per_job - 1;

disp(['Running simulations ',num2str(start_end(1)), ' to ', num2str(start_end(2))]);

%Loop through points
for i = start_end(1):start_end(2)
    
    
    %Only run simulations that are required
    if i <= size(par_table,1)
        
        %Read in default parameter structure
        load('parameters/default_parameter_struct.mat');
        
        %Replace parameters
        non_id_vars = par_table.Properties.VariableNames(1:(end-1));
        matching_row = par_table.id == i;
        for j = 1:length(non_id_vars)
            par.(non_id_vars{j}) = par_table{matching_row,non_id_vars{j}};
        end

        %If multiple Keq or K exist, convert
        if isfield(par,'Keq_1')
            par.Keq = [par.Keq_1,par.Keq_2];
        end
        if isfield(par,'K_1')
            par.K = [par.K_1,par.K_2];
        end
        if isfield(par,'alpha_1')
            par.K = [par.alpha_1,par.alpha_2];
        end

        true_ind = i - start_end(1) + 1;
        disp(['Running simulation ',num2str(i),'.'])
        
        %Run simulation
        [ET_cell{true_ind},~,growth_cell{true_ind},~,~,simu_cell{true_ind},num_auto_peak_cell{true_ind}]...
            = chemostat_invasion_ode15s();
    end
    
end

%Generate results file
results_file = [results_folder,'/','sim_id_',num2str(batch_id),'.mat'];

save(results_file,'ET_cell','growth_cell','simu_cell','start_end','num_auto_peak_cell');

disp(['Raw results saved to ',results_file,'.'])


end
