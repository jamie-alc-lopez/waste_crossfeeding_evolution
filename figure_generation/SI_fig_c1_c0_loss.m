%% This script generates the SI figure plotting the loss of c1 in the P1C 
% and P1P2 consortia. These values are generated from numerical simulation

%First, load and process data

clear;clc

%Laod in numerical data and default simulation parameters
load('fig_data/P1C_vs_P1P2_linear_sweep_simulation_results.mat')
load('../model_code/parameters/default_parameter_struct.mat')

par_tables{1} = parameter_table(parameter_table.Keq_1 == 1,:);
par_tables{2} = parameter_table(parameter_table.Keq_1 == 10,:);

%Reference enzyme compositions
P1_enz = [1 0 1 1 0];
P2_enz = [1 1 1 0 1];
C_enz = [0 1 0 1 1];
min_E = 1e-6;

%Reformat data
for k = 1:length(par_tables)

    %Identify consortia composition
    only_P1P2_ind{k} = par_tables{k}.only_P1P2 == 1;
    only_P1C_ind{k} = par_tables{k}.only_P1C == 1;
    consortium_ind{k} = cellfun(@(x) size(x,1),par_tables{k}.strains);
    consortium_ind{k} = consortium_ind{k} == 2;

    %Loop through table
    for i = 1:size(par_tables{k},1)

        %Extract the min c0 and c1/c2 values
        par_tables{k}.min_c0(i) = par_tables{k}.simu{i}.c(end,1);
        par_tables{k}.c1(i) = par_tables{k}.simu{i}.c(end,2);
        par_tables{k}.c2(i) = par_tables{k}.simu{i}.c(end,3);

        %Get strain abundances and 
        if consortium_ind{k}(i)
            %Locate strain types
            [~,par_tables{k}.P1_ind{i}] = intersect(par_tables{k}.strains{i}>min_E,P1_enz,'rows');
            [~,par_tables{k}.P2_ind{i}] = intersect(par_tables{k}.strains{i}>min_E,P2_enz,'rows');
            [~,par_tables{k}.C_ind{i}] = intersect(par_tables{k}.strains{i}>min_E,C_enz,'rows');

            %Extract populations and internal concentrations
            %P1
            par_tables{k}.P1_rho(i) = par_tables{k}.simu{i}.rho(par_tables{k}.P1_ind{i});
            par_tables{k}.P1_c0(i) = par_tables{k}.simu{i}.c(par_tables{k}.P1_ind{i},1);
            par_tables{k}.P1_c1(i) = par_tables{k}.simu{i}.c(par_tables{k}.P1_ind{i},2);
            par_tables{k}.P1_c2(i) = par_tables{k}.simu{i}.c(par_tables{k}.P1_ind{i},3);

            %C
            if only_P1C_ind{k}(i)
                par_tables{k}.C_rho(i) = par_tables{k}.simu{i}.rho(par_tables{k}.C_ind{i});
                par_tables{k}.C_c0(i) = par_tables{k}.simu{i}.c(par_tables{k}.C_ind{i},1);
                par_tables{k}.C_c1(i) = par_tables{k}.simu{i}.c(par_tables{k}.C_ind{i},2);
                par_tables{k}.C_c2(i) = par_tables{k}.simu{i}.c(par_tables{k}.C_ind{i},3);

            %P2
            elseif only_P1P2_ind{k}(i)
                par_tables{k}.P2_rho(i) = par_tables{k}.simu{i}.rho(par_tables{k}.P2_ind{i});
                par_tables{k}.P2_c0(i) = par_tables{k}.simu{i}.c(par_tables{k}.P2_ind{i},1);
                par_tables{k}.P2_c1(i) = par_tables{k}.simu{i}.c(par_tables{k}.P2_ind{i},2);
                par_tables{k}.P2_c2(i) = par_tables{k}.simu{i}.c(par_tables{k}.P2_ind{i},3);
            end
        end

    end

    %Process into reference vectors for plotting
    supply_vec{k} = par_tables{k}.supply(only_P1P2_ind{k});
    valid_ratio_ind{k} = consortium_ind{k}(only_P1C_ind{k}) & consortium_ind{k}(only_P1P2_ind{k});
    valid_supply_vec{k} = supply_vec{k}(valid_ratio_ind{k});

end



%% Plot the c1 loss rates


%Generate empty figure
newfigure(3.42/2, (1.5/3)*3.42/2);
gap = [0.2,0.2]; %height, width
marg_h = [0.2,0.1]; %lower upper
marg_w = [0.15,0.05]; %left right
FontSize = 4;
LineWidth = 0.4;
Figax = tight_subplot(1,2,gap,marg_h,marg_w);
colors = [220 38 127;...
    220 38 127]./255;
colors(2,:) = colors(2,:)*0.6;
dashtype = {'-',':','-',':'};

%Loop through Keq values and plot


for k = 1:2

    %Separate by P1C or P1P2
    table_k = par_tables{k};
    table_k_P1C = table_k(table_k.only_P1C==1,:);
    table_k_P1P2 = table_k(table_k.only_P1P2==1,:);

    %Compute losses of c1
    P1C_c1_loss = table_k_P1C.c1.*table_k_P1C.D + ... %Extracellular dilution (c1e*D)
        table_k_P1C.P1_c1.*par.VcVC.*table_k_P1C.P1_rho.*table_k_P1C.D + ... %Dilution of P1 (c1_P1*rV*rho_P1*D)
        table_k_P1C.C_c1.*par.VcVC.*table_k_P1C.C_rho.*table_k_P1C.D; %Dilution of C (c1_C*rV*rho_C*D)

    P1P2_c1_loss = table_k_P1P2.c1.*table_k_P1P2.D + ... %Extracellular dilution (c1e*D)
        table_k_P1P2.P1_c1.*par.VcVC.*table_k_P1P2.P1_rho.*table_k_P1P2.D + ... %Dilution of P1 (c1_P1*rV*rho_P1*D)
        table_k_P1P2.P2_c1.*par.VcVC.*table_k_P1P2.P2_rho.*table_k_P1P2.D; %Dilution of P2 (c1_P2*rV*rho_P2*D)

    %Compute losses of c0
    P1C_c0_loss = table_k_P1C.min_c0.*table_k_P1C.D + ... %Extracellular dilution (c0e*D)
        table_k_P1C.P1_c0.*par.VcVC.*table_k_P1C.P1_rho.*table_k_P1C.D; %Dilution of P1 (c0_P1*rV*rho_P1*D)

    P1P2_c0_loss = table_k_P1P2.min_c0.*table_k_P1P2.D + ... %Extracellular dilution (c0e*D)
        table_k_P1P2.P1_c0.*par.VcVC.*table_k_P1P2.P1_rho.*table_k_P1P2.D + ... %Dilution of P1 (c0_P1*rV*rho_P1*D)
        table_k_P1P2.P2_c0.*par.VcVC.*table_k_P1P2.P2_rho.*table_k_P1P2.D; %Dilution of P2 (c0_P2*rV*rho_P2*D)


    ind = valid_ratio_ind{k};

    %Plot
    axes(Figax(1))
    hold on
    plot(valid_supply_vec{k},(P1C_c1_loss(ind)-P1P2_c1_loss(ind))./valid_supply_vec{k},...
        dashtype{1},'LineWidth',LineWidth,'Color',colors(k,:))

    axes(Figax(2))
    hold on
    plot(valid_supply_vec{k},(P1C_c0_loss(ind)-P1P2_c0_loss(ind))./valid_supply_vec{k},...
        dashtype{1},'LineWidth',LineWidth,'Color',colors(k,:))


end


ylabel1_cell = {'Norm. difference in $c_1$ loss','Norm. difference in $c_0$ loss'};
label_cell = {'\textbf{A}','\textbf{B}'};

for i = 1:2
    axes(Figax(i))

    %Set up axes
    xlim([1e-8,1e-5])
    ylim([0,0.15])
    xticks([1e-8,1e-5]);
    yticks([0,0.05,0.1,0.15]);
    xticklabels({'$10^{-8}$','$10^{-5}$'})
    yticklabels({'0','0.05','0.1','0.15'})

    set(gca,'XScale','log')
    set(gca,'XMinorTick','Off')
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',FontSize)

    ylabel({ylabel1_cell{i},'rate (P1+C minus P1+P2)'},'Interpreter','latex','FontSize',FontSize)
    xlabel('Supply rate, $s_0$','Interpreter','latex','FontSize',FontSize)

    text(-0.42,1.09,label_cell{i},'Interpreter','latex','Units','normalized','FontSize',4)
    set(gca,'Clipping','off')

end

%Make legend
axes(Figax(2))
leg_top = 0.14;
leg_left = 1e-7;
mylines = {'-','-'};
leg_colors = {colors(2,:),colors(1,:)};
labels = {'$K_\textrm{eq} = 10$','$K_\textrm{eq} = 1$'};
line_length = 3e-7;
spacing_x = 1e-7;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',leg_colors{k})
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.016;
end



print(gcf, '-dpng','figs/SI_fig_c1_c0_loss.png','-r1200','-painters');
