%% This script generates the parameter tables for the HPC simulation runs
% used in the main text

%% Non-saturating kinetics with osmolarity, Keq 1

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[1,1]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 1;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/osm_sweep_Keq_1_no_sat.csv','Delimiter',',')


%% Non-saturating kinetics with osmolarity, Keq 10

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[10,10]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 1;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/osm_sweep_Keq_10_no_sat.csv','Delimiter',',')



%% Saturating kinetics with no osmolarity, Keq 1

clear;clc

n_ind = 141;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 0;
        parameter_table{tn,'fast_pseudo'} = 0;
        parameter_table{tn,'Keq'} = {[1,1]};
        parameter_table{tn,'K'} = {[1,1]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 10;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq_1_sat.csv','Delimiter',',')

%% Saturating kinetics with no osmolarity, Keq 10

clear;clc

n_ind = 141;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 0;
        parameter_table{tn,'fast_pseudo'} = 0;
        parameter_table{tn,'Keq'} = {[10,10]};
        parameter_table{tn,'K'} = {[1,1]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 10;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq_10_sat.csv','Delimiter',',')


%% Linear kinetics with no osmolarity, Keq = 1

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[1,1]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq_1_no_sat.csv','Delimiter',',')

%% Linear kinetics with no osmolarity, Keq = 10

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[10,10]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq_10_no_sat.csv','Delimiter',',')


%% Linear kinetics, no osmolarity, variable Keq

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[1,0.1]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq1_1_Keq2_0p1_no_sat.csv','Delimiter',',')


%% Linear kinetics, no osmolarity, Keq 1, alpha_2 = 0

clear;clc

n_ind = 200;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo'};
n_par = length(table_columns);
n_tot = n_ind*n_ind;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_range = logspace(-6,-3,n_ind);
s_range = logspace(-8,-5,n_ind);

tn = 1;
for i = 1:n_ind
    for j = 1:n_ind
        parameter_table{tn,'D'} = D_range(i);
        parameter_table{tn,'supply'} = s_range(j);
        parameter_table{tn,'no_saturation'} = 1;
        parameter_table{tn,'fast_pseudo'} = 1;
        parameter_table{tn,'Keq'} = {[1,1]};
        parameter_table{tn,'alpha'} = {[1.82e-1,0]};
        tn = tn + 1;
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 50;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,'hpc_tables/no_osm_sweep_Keq_1_alpha2_0_no_sat.csv','Delimiter',',')


%% Small sweep to compare optimal P1/P2 and P1C

clear;clc

n_ind = 200;
n_restrictions = 3;
n_Keq = 2;

table_columns = {'D','supply','osm_coeff','no_saturation','fast_pseudo',...
    'only_P1P2','only_P1C'};
n_par = length(table_columns);
n_tot = n_ind*n_restrictions*n_Keq;

parameter_table = nan(n_tot,n_par);

parameter_table = array2table(parameter_table,'VariableNames',table_columns);

D_val = 1e-6;
s_range = logspace(-8,-5,n_ind);
Keq_range = [1,10];


tn = 1;
for w = 1:n_Keq
    for i = 1:n_restrictions
        for j = 1:n_ind

            parameter_table{tn,'only_P1P2'} = 0;
            parameter_table{tn,'only_P1C'} = 0;
            if i == 1
                parameter_table{tn,'only_P1P2'} = 1;
            elseif i == 2
                parameter_table{tn,'only_P1C'} = 1;
            end

            parameter_table{tn,'D'} = D_val;
            parameter_table{tn,'supply'} = s_range(j);
            parameter_table{tn,'no_saturation'} = 1;
            parameter_table{tn,'fast_pseudo'} = 1;
            parameter_table{tn,'Keq'} = {[Keq_range(w),Keq_range(w)]};
            tn = tn + 1;

        end
    end
end

parameter_table{:,'id'} = (1:n_tot)';
parameter_table.osm_coeff(:) = 0;

run_per_job = 10;

disp(['These ', num2str(size(parameter_table,1)),' simulations will require ', num2str(ceil(size(parameter_table,1)/run_per_job)),' jobs at ', num2str(run_per_job),' runs per job.']); 

writetable(parameter_table,['hpc_tables/P1C_vs_P1P2_linear_sweep.csv'],'Delimiter',',')
