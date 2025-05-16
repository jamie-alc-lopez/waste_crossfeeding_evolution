%% This script generates figure 3 

clear;clc
close all

%Laod in numerical data
load('fig_data/P1C_vs_P1P2_linear_sweep_simulation_results.mat')

par_tables{1} = parameter_table(parameter_table.Keq_1 == 1,:);
par_tables{2} = parameter_table(parameter_table.Keq_1 == 10,:);

%Reformat data
for k = 1:length(par_tables)

    for i = 1:size(par_tables{k},1)
        par_tables{k}.min_c0(i) = par_tables{k}.simu{i}.c(end,1);
    end
    only_P1P2_ind{k} = par_tables{k}.only_P1P2 == 1;
    only_P1C_ind{k} = par_tables{k}.only_P1C == 1;
    consortium_ind{k} = cellfun(@(x) size(x,1),par_tables{k}.strains);
    consortium_ind{k} = consortium_ind{k} == 2;
    supply_vec{k} = par_tables{k}.supply(only_P1P2_ind{k});
    c0e_ratio{k} =par_tables{k}.min_c0(only_P1C_ind{k})./par_tables{k}.min_c0(only_P1P2_ind{k});
    valid_ratio_ind{k} = consortium_ind{k}(only_P1C_ind{k}) & consortium_ind{k}(only_P1P2_ind{k});

    valid_supply_vec{k} = supply_vec{k}(valid_ratio_ind{k});
    valid_c0e_ratio{k} = c0e_ratio{k}(valid_ratio_ind{k});

end

%Initialize figure
newfigure(3.42/2, (1.5/3)*3.42/2);
gap = [0.2,0.2]; %height, width
marg_h = [0.2,0.1]; %lower upper
marg_w = [0.15,0.05]; %left right
FontSize = 4;
LineWidth = 0.4;
Figax = tight_subplot(1,2,gap,marg_h,marg_w);

colors = repmat([220 38 127],2,1)./255;
colors(2,:) = colors(2,:)*0.6;

%Plot panel B
axes(Figax(1))
hold on
dashes = {'-','-'};
plot([supply_vec{1}(1),supply_vec{1}(end)],[1,1],':','LineWidth',LineWidth,'color',0.5*[1,1,1])

for k = 1:length(valid_supply_vec)
    plot(valid_supply_vec{k},valid_c0e_ratio{k},...
        dashes{k},'Color',colors(k,:),'LineWidth',LineWidth)
end

xlim([supply_vec{1}(1),supply_vec{1}(end)])
ylim([0.95,1.3])
xticks([supply_vec{1}(1),supply_vec{1}(end)]);
yticks([1,1.1,1.2,1.3]);
xticklabels({'$10^{-8}$','$10^{-5}$'})
yticklabels({'1','1.1','1.2','1.3'})

set(gca,'XScale','log')
set(gca,'XMinorTick','Off')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)

ylabel({'Minimum P1+C to','P1+P2 $c_{0,e}$ ratio'},'Interpreter','latex','FontSize',FontSize)
xlabel('Supply rate, $s_0$','Interpreter','latex','FontSize',FontSize)

leg_top = 1.3;
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
    curr_y = curr_y-0.04;
end
text(-0.42,1.06,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'Clipping','off')


%% Plot internal concentrations of the consortia (panel C)

%Loop through results and extract the relevant concentration values
P1_enz = [1 0 1 1 0];
P2_enz = [1 1 1 0 1];
C_enz = [0 1 0 1 1];
min_E = 1e-6;
for i = 1:size(par_tables{1},1)
    if consortium_ind{1}(i)

        %Locate strain types
        [~,par_tables{1}.P1_ind{i}] = intersect(par_tables{1}.strains{i}>min_E,P1_enz,'rows');
        [~,par_tables{1}.P2_ind{i}] = intersect(par_tables{1}.strains{i}>min_E,P2_enz,'rows');
        [~,par_tables{1}.C_ind{i}] = intersect(par_tables{1}.strains{i}>min_E,C_enz,'rows');

        %Extract internal concentrations
        par_tables{1}.int_c0(i) = par_tables{1}.simu{i}.c(par_tables{1}.P1_ind{i},1);
        if only_P1C_ind{1}(i)
            par_tables{1}.int_c2(i) = par_tables{1}.simu{i}.c(par_tables{1}.C_ind{i},3);
        elseif only_P1P2_ind{1}(i)
            par_tables{1}.int_c2(i) = par_tables{1}.simu{i}.c(par_tables{1}.P2_ind{i},3);
        end

    end
end

%Process into plotting vectors
c0_P1P2 = par_tables{1}.int_c0(only_P1P2_ind{1});
c0_P1P2 = c0_P1P2(valid_ratio_ind{1});

c0_P1C = par_tables{1}.int_c0(only_P1C_ind{1});
c0_P1C = c0_P1C(valid_ratio_ind{1});

c2_P1P2 = par_tables{1}.int_c2(only_P1P2_ind{1});
c2_P1P2 = c2_P1P2(valid_ratio_ind{1});

c2_P1C = par_tables{1}.int_c2(only_P1C_ind{1});
c2_P1C = c2_P1C(valid_ratio_ind{1});

xvec = valid_supply_vec{1};
%norm_vec = ones(size(xvec)); 
norm_vec = xvec./par_tables{1}.D(1);


%Plot the figure
axes(Figax(2))

colors1 = [64 176 166;64 176 166]/255;
colors1 = colors1.*[0.9;1.4];
dashtype = {'-',':','-',':'};
sim_ind = 1;

hold on

plot(xvec,c0_P1P2./norm_vec,dashtype{1},'LineWidth',LineWidth,'Color',colors1(1,:))
plot(xvec,c2_P1P2./norm_vec,dashtype{2},'LineWidth',LineWidth,'Color',colors1(1,:))
plot(xvec,c0_P1C./norm_vec,dashtype{3},'LineWidth',LineWidth,'Color',colors1(2,:))
plot(xvec,c2_P1C./norm_vec,dashtype{4},'LineWidth',LineWidth,'Color',colors1(2,:))

y_lim_vec = [0,0.7];
ylim(y_lim_vec)
yticks(y_lim_vec);
yticklabels({'0','0.7'})

xlim([supply_vec{1}(1),supply_vec{1}(end)])
xticks([supply_vec{1}(1),supply_vec{1}(end)]);
xticklabels({'$10^{-8}$','$10^{-5}$'})

set(gca,'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)

set(gca,'XScale','log')
%set(gca,'YScale','log')
set(gca,'XMinorTick','Off')

ylabel({'Normalized intracellular','concentration, $(c_i^* D)/s_0$'},'Interpreter','latex','FontSize',FontSize)
xlabel('Supply rate, $s_0$','Interpreter','latex','FontSize',FontSize)

leg_top = y_lim_vec(end);
leg_left = 1e-7;
leg_colors = {colors1(1,:),colors1(1,:),colors1(2,:),colors1(2,:)};
labels = {'$c_0^*$ (\textbf{P1}+P2)','$c_2^*$ (P1+\textbf{P2})',...
    '$c_0^*$ (\textbf{P1}+C)','$c_2^*$ (P1+\textbf{C})'};
line_length = 0.3e-6;
spacing_x = 1e-7;
curr_y = leg_top;
order = [3 1 4 2];
for k = 1:length(labels)
    w = order(k);
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],dashtype{w},'Color',leg_colors{w})
    text(leg_left + line_length + spacing_x,curr_y,labels{w},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.075;
end

text(-0.45,1.06,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)



print(gcf, '-dpng','figs/fig3BC.png','-r1200','-painters');
