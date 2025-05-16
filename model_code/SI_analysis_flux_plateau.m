clear;clc

%This script examines the growth of P1 strains near the second transition
%in the cross-feeding phase

sim_file = '../figure_generation/fig_data/osm_sweep_Keq_1_no_sat_simulation_results.mat';

[lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file);
D_range = [sr.D_range(1),sr.D_range(end)];
S_range = [sr.supply_range(1),sr.supply_range(end)];

ifn = @(x) x(:,22); %57 for old fig, 22 for new fig

extra_c1 = ifn(sr.extra_c_cell{2});
P1_E1 = ifn(sr.P1_E_cell{1});

ET_cell = ifn(sr.raw_E_cell);
c_cell = ifn(sr.raw_c_cell);
simu_cell = ifn(sr.simu_cell);
D_vec = ifn(sr.D_mat);
supply_vec = ifn(sr.supply_mat);

%Pre-transition, during transition, and post-transition comparisons
comparison_pairs = [161 164; 164 167; 167 170];

E_cutoff = 1e-6;
strain_types = [1,0,1,1,0];

VariableNames = {'dilution1','supply_rate1','dilution2','supply_rate2',...
    'g_P1_1_in_env1','g_P1_1_in_env2','g_P1_2_in_env1','g_P1_2_in_env2',...
    'osm_P1_1_in_env1','osm_P1_1_in_env2','osm_P1_2_in_env1','osm_P1_2_in_env2',...
    'flux_P1_1_in_env1','flux_P1_1_in_env2','flux_P1_2_in_env1','flux_P1_2_in_env2'};
RowNames = {'pre-to-pre','pre-to-post','post-to-post'};
out_table = array2table(nan(size(comparison_pairs,1),length(VariableNames)),...
    'VariableNames',VariableNames,'RowNames',RowNames);

%Loop through comparison pairs
for i = 1:size(comparison_pairs,1)

    ind1 = comparison_pairs(i,1);
    ind2 = comparison_pairs(i,2);

    out_table.dilution1(i) = D_vec(ind1);
    out_table.dilution2(i) = D_vec(ind2);

    out_table.supply_rate1(i) = supply_vec(ind1);
    out_table.supply_rate2(i) = supply_vec(ind2);

    ET1 = ET_cell{ind1};
    ET2 = ET_cell{ind2};

    [~,temp_ind]=intersect(ET1>E_cutoff,strain_types,'rows');
    P1_1 = ET1(temp_ind,:);

    [~,temp_ind]=intersect(ET1>E_cutoff,strain_types,'rows');
    P1_2 = ET2(temp_ind,:);

    simu1 = simu_cell{ind1};
    simu2 = simu_cell{ind2};

    [out_table.g_P1_1_in_env1(i),~,~,...
        out_table.flux_P1_1_in_env1(i),out_table.osm_P1_1_in_env1(i)] = ...
        compute_pseudo_SS_ode15s(P1_1,simu1);

    [out_table.g_P1_1_in_env2(i),~,~,...
        out_table.flux_P1_1_in_env2(i),out_table.osm_P1_1_in_env2(i)] = ...
        compute_pseudo_SS_ode15s(P1_1,simu2);

    [out_table.g_P1_2_in_env1(i),~,~,...
        out_table.flux_P1_2_in_env1(i),out_table.osm_P1_2_in_env1(i)] = ...
        compute_pseudo_SS_ode15s(P1_2,simu1);

    [out_table.g_P1_2_in_env2(i),~,~,...
        out_table.flux_P1_2_in_env2(i),out_table.osm_P1_2_in_env2(i)] = ...
        compute_pseudo_SS_ode15s(P1_2,simu2);


end

out_table.Variables = round(out_table.Variables,3,'significant');

%writetable(out_table,'second_transition_growth_analysis.xlsx',...
%    'WriteRowNames',true,'WriteVariableNames',true);