%% This script generates the SI figure that shows a test for hysteresis

clear;clc

global par

%Perform hysteresis test on subset of environments from larger set of
%simulations
load('../model_code/parameters/default_parameter_struct.mat');
table_file = '../model_code/parameters/hpc_tables/osm_sweep_Keq_1_no_sat.csv';

environments = [131 132 133 134 135 134 133 132 131];
par.extol = 1e3;
[ET_store,simu_store] = test_for_hysteresis(environments,table_file);

save('fig_data/hysteresis_test.mat')

%% Analyze hysteresis test

clear;clc

%Load hysteresis data and full simulation data
load('fig_data/hysteresis_test.mat')
load('fig_data/osm_sweep_Keq_1_no_sat_simulation_results.mat')

strain_types = [1 0 1 1 0; 1 1 1 0 1; 0 1 0 1 1; 1 1 1 1 1];
phases = {1,2,[1;2],[1;3]};
for i = 1:length(environments)
    
    current_ET = ET_store{i};
    simplified_ET = unique(current_ET > E_cutoff,'rows');

    environment_ind = parameter_table.id == environments(i);

    if enzyme_only
        [~,types_present] = intersect(strain_types(1:end-1,1:2),simplified_ET(:,1:2),'rows');
    else
        [~,types_present] = intersect(strain_types,simplified_ET,'rows');
    end
    types_present = sort(types_present,'ascend');
    expected_strain_types{i} = parameter_table.strain_types{environment_ind};
    present_strain_types{i} = types_present;
    
    test_phase_matrix(1,i) = find(cellfun(@(x) isequal(x,expected_strain_types{i}),phases));
        test_phase_matrix(2,i) = find(cellfun(@(x) isequal(x,types_present),phases));

    test_supply_vec(i) = parameter_table.supply(environment_ind);
    test_D_vec(i) = parameter_table.D(environment_ind);

    
end


figure
imagesc(test_phase_matrix,'AlphaData',~isnan(test_phase_matrix),[0,4])
yticks([1,2])
xticklabels(arrayfun(@(x) num2str(round(x*1e7,2)), test_supply_vec,'UniformOutput',false))
yticklabels({'Abiotic','Stepwise'})

xlabel('Supply rate, $s_0$ ($\times 10^{-7}$)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)

pos = get(gcf,'Position');
pos(4) = 0.5*pos(4);
set(gcf,'Position',pos);

full_phase = [expected_strain_types;present_strain_types];
hold on
for i = 1:size(test_phase_matrix,1)
    plot([9.5,0.5],[i+0.5,i+0.5],'k-')
    for j = 1:size(test_phase_matrix,2)
        plot([j+0.5,j+0.5],[3.5,0.5],'k-')
        if isequal(full_phase{i,j},[1;3])
            text(j,i,'P1+C','HorizontalAlignment','center','Interpreter','latex')
        else
            text(j,i,'P1+P2','HorizontalAlignment','center','Interpreter','latex')
        end
    end
end

print(gcf, '-dpng',['figs/SI_fig_hysteresis_test.png'],'-r700');
