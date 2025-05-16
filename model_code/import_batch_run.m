function [parameter_table,sr] = import_batch_run(sim_name,num_jobs,run_per_job,E_cutoff)

%This script imports the results of a large-scale simulation run

%Import parameter table
parameter_table = readtable(['../../parameters/hpc_tables/',sim_name,'.csv']);

%Get default parameter struct
load('../../parameters/default_parameter_struct.mat');

%Define strain types
strain_types = [1 0 1 1 0; 1 1 1 0 1; 0 1 0 1 1; 1 1 1 1 1];

%Define phases of interest
phases = {1,2,[1;2],[1;3],3,[2;3]};

%Initialized result struct
sr = struct();

%Extract parameter ranges
sr.supply_range = unique(parameter_table.supply);
sr.D_range = unique(parameter_table.D);

%Generate storage structures
sr.phase_cell = cell(length(sr.supply_range),length(sr.D_range));
sr.simu_cell = sr.phase_cell;
sr.raw_c_cell = sr.phase_cell;
sr.raw_E_cell = sr.phase_cell;
sr.raw_rho_cell = sr.phase_cell;
sr.phase_matrix = nan(length(sr.supply_range),length(sr.D_range));
sr.supply_mat = zeros(size(sr.phase_matrix));
sr.D_mat = sr.supply_mat;
sr.degradation_mat = sr.supply_mat;
sr.V_mat = sr.supply_mat;
sr.C_mat = sr.supply_mat;
sr.P1_mat = sr.supply_mat;
sr.P2_mat = sr.supply_mat;
sr.C_flux_mat = sr.supply_mat;
sr.P1_flux_mat = sr.supply_mat;
sr.P2_flux_mat = sr.supply_mat;

%Make concentration storage cell
c_cell = cell(3,1);
for i = 1:3
    c_cell{i} = sr.supply_mat;
end
sr.V_c_cell = c_cell;
sr.C_c_cell = c_cell;
sr.P1_c_cell = c_cell;
sr.P2_c_cell = c_cell;
sr.extra_c_cell = c_cell;
sr.burden_mat = sr.supply_mat;

%Make enzyme storage cell
E_cell = cell(5,1);
for i = 1:5
    E_cell{i} = sr.supply_mat;
end
sr.E_cell = E_cell;
sr.C_E_cell = E_cell;
sr.P1_E_cell = E_cell;
sr.P2_E_cell = E_cell;

%Loop through files and fill storage structures
for i = 1:num_jobs
    
    %Load file if it exists
    sim_file = ['../results/',sim_name,'/sim_id_',num2str(i),'.mat'];
    if isfile(sim_file)
        load(sim_file)
        
        for j = 1:length(simu_cell)
            true_ind = run_per_job*(i-1) + j;
            current_ET = ET_cell{j};
            simu = simu_cell{j};
            opt_growth = growth_cell{j};
            
            current_ET = current_ET(simu.rho > 1000,:);
            
            %Identify which types of strains exist at steady state
            parameter_table.strains{true_ind} = current_ET;
            simplified_ET = unique(current_ET > E_cutoff,'rows');
            [~,types_present] = intersect(strain_types,simplified_ET,'rows');
            types_present = sort(types_present,'ascend');
            parameter_table.strain_types{true_ind} = types_present;
            
            %Determine to indices of this simulation's parameters
            [~,ind1] = intersect(sr.supply_range,parameter_table.supply(true_ind));
            [~,ind2] = intersect(sr.D_range,parameter_table.D(true_ind));
            sr.phase_cell{ind1,ind2} = parameter_table.strain_types{true_ind};
            sr.supply_mat(ind1,ind2) = parameter_table.supply(true_ind);
            sr.D_mat(ind1,ind2) = parameter_table.D(true_ind);
            
            %Find phase matching list of strains
            if isempty(types_present)
                sr.phase_matrix(ind1,ind2) = 0;
            elseif ~isempty(num_auto_peak_cell{j})
                sr.phase_matrix(ind1,ind2) = -1;
            else
                matching_phase = find(cellfun(@(x) isequal(x,types_present),phases));
                if ~isempty(matching_phase)
                    sr.phase_matrix(ind1,ind2) = matching_phase;
                else
                    sr.phase_matrix(ind1,ind2) = -2;
                end
            end
            parameter_table.phase_index(true_ind) = sr.phase_matrix(ind1,ind2);
            parameter_table.batch_id(true_ind) = i; 
            parameter_table.ind1(true_ind) = ind1;
            parameter_table.ind2(true_ind) = ind2;

            %Record end simulation state
            sr.raw_E_cell{ind1,ind2} = current_ET;
            sr.raw_c_cell{ind1,ind2} = simu.c;
            sr.extra_c_cell{1}(ind1,ind2) = simu.c(end,1);
            sr.extra_c_cell{2}(ind1,ind2) = simu.c(end,2);
            sr.extra_c_cell{3}(ind1,ind2) = simu.c(end,3);
            parameter_table.simu{true_ind} = simu;
            
            %Store analysis results
            sr.raw_rho_cell{ind1,ind2} = simu.rho;
            
            %Store simu object
            sr.simu_cell{ind1,ind2} = simu;

            %Store strain-specific data
            for k = 1:size(current_ET,1)
                simplified_ET = current_ET(k,:) > E_cutoff;
                [~,species_type] = intersect(strain_types,simplified_ET,'rows');
                if species_type == 1
                    sr.P1_mat(ind1,ind2) = sr.P1_mat(ind1,ind2) + simu.rho(k);
                    sr.P1_c_cell{1}(ind1,ind2) = simu.c(k,1);
                    sr.P1_c_cell{2}(ind1,ind2) = simu.c(k,2);
                    sr.P1_c_cell{3}(ind1,ind2) = simu.c(k,3);
                    sr.P1_E_cell{1}(ind1,ind2) = current_ET(k,1);
                    sr.P1_E_cell{2}(ind1,ind2) = current_ET(k,2);
                    sr.P1_E_cell{3}(ind1,ind2) = current_ET(k,3);
                    sr.P1_E_cell{4}(ind1,ind2) = current_ET(k,4);
                    sr.P1_E_cell{5}(ind1,ind2) = current_ET(k,5);
                    sr.P1_flux_mat(ind1,ind2) = par.VcVC*par.beta*current_ET(k,3)*(simu.c(end,1) - simu.c(k,1));

                elseif species_type == 2
                    sr.P2_mat(ind1,ind2) = sr.P2_mat(ind1,ind2) + simu.rho(k);
                    sr.P2_c_cell{1}(ind1,ind2) = simu.c(k,1);
                    sr.P2_c_cell{2}(ind1,ind2) = simu.c(k,2);
                    sr.P2_c_cell{3}(ind1,ind2) = simu.c(k,3);
                    sr.P2_E_cell{1}(ind1,ind2) = current_ET(k,1);
                    sr.P2_E_cell{2}(ind1,ind2) = current_ET(k,2);
                    sr.P2_E_cell{3}(ind1,ind2) = current_ET(k,3);
                    sr.P2_E_cell{4}(ind1,ind2) = current_ET(k,4);
                    sr.P2_E_cell{5}(ind1,ind2) = current_ET(k,5);
                    sr.P2_flux_mat(ind1,ind2) = par.VcVC*par.beta*current_ET(k,3)*(simu.c(end,1) - simu.c(k,1));
                elseif species_type == 3
                    sr.C_mat(ind1,ind2) = sr.C_mat(ind1,ind2) + simu.rho(k);
                    sr.C_c_cell{1}(ind1,ind2) = simu.c(k,1);
                    sr.C_c_cell{2}(ind1,ind2) = simu.c(k,2);
                    sr.C_c_cell{3}(ind1,ind2) = simu.c(k,3);                    
                    sr.C_E_cell{1}(ind1,ind2) = current_ET(k,1);
                    sr.C_E_cell{2}(ind1,ind2) = current_ET(k,2);
                    sr.C_E_cell{3}(ind1,ind2) = current_ET(k,3);
                    sr.C_E_cell{4}(ind1,ind2) = current_ET(k,4);
                    sr.C_E_cell{5}(ind1,ind2) = current_ET(k,5);
                    sr.C_flux_mat(ind1,ind2) = par.VcVC*par.beta*current_ET(k,4)*(simu.c(end,2) - simu.c(k,2));
                elseif species_type == 4
                    sr.V_mat(ind1,ind2) = sr.V_mat(ind1,ind2) + simu.rho(k);
                    sr.V_c_cell{1}(ind1,ind2) = simu.c(k,1);
                    sr.V_c_cell{2}(ind1,ind2) = simu.c(k,2);
                    sr.V_c_cell{3}(ind1,ind2) = simu.c(k,3);
                    sr.V_E_cell{1}(ind1,ind2) = current_ET(k,1);
                    sr.V_E_cell{2}(ind1,ind2) = current_ET(k,2);
                    sr.V_E_cell{3}(ind1,ind2) = current_ET(k,3);
                    sr.V_E_cell{4}(ind1,ind2) = current_ET(k,4);
                    sr.V_E_cell{5}(ind1,ind2) = current_ET(k,5);
                    disp('Versatilist strain detected') 
                else
                    disp('Unclassified strain detected')
                end
            end
            
            %Compute mean burdens
            rho = simu.rho;
            c = simu.c;
            burden = 0;
            for k = 1:length(rho)
                burden = burden + rho(k)*sum(c(k,:));
            end
            sr.burden_mat(ind1,ind2) = burden/sum(rho);

            
            
        end
    else
        disp(['File for simulation batch ',num2str(i),' not found.'])
    end
end

%Compute degradation fluxes
%sr.degradation_mat = sr.D_mat.*cell2mat(cellfun(@(x) sum(x(end,:)),sr.raw_c_cell,'UniformOutput',false));


end
