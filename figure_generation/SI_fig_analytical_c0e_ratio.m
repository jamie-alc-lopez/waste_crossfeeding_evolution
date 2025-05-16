%% This script generates the SI figure plotting the P1C to P1P2 c0e ratio 
% in the model neglecting growth dilution

% This first part computes the c0e values, data is pregenerated it can be
% skipped

clear;clc

Vc = 2e-15;
alpha= 1.8e-1;
beta = 2e2;
Keq_vals = [1,10];
D_vals = 1e-6;
gm = 5e-4;
n_points_S = 80;
S_vec = logspace(-8,-5,n_points_S);

for i = 1:length(Keq_vals)
    for j = 1:length(D_vals)
        for k = 1:length(S_vec)
            [ESS_ind_mat(i,j,k),ESS_type_cell{i,j,k},...
                optimal_ET_cell{i,j,k},optimal_c0e_cell{i,j,k},ss_val_cell{i,j,k},all_c0e{i,j,k},xopt_cell{i,j,k}]...
                = compute_ESS_consortium(S_vec(k),alpha,alpha,beta,D_vals(j),Vc,Keq_vals(i),Keq_vals(i),gm);
            P1C_P1P2_ratio_analytical(i,j,k) = P1C_P1P2_c0e_ratio(S_vec(k),alpha,alpha,beta,D_vals(j),Keq_vals(i),Keq_vals(i),gm);
        end
    end
end

save('fig_data/analytical_P1C_to_P1P2_comparison.mat')

%% Plot SI fig panel C - internal concentrations
clear;clc

load('fig_data/analytical_P1C_to_P1P2_comparison.mat')

newfigure(3.42/2, (1.5/3)*3.42/2);
gap = [0.2,0.2]; %height, width
marg_h = [0.2,0.1]; %lower upper
marg_w = [0.15,0.05]; %left right
FontSize = 4;
LineWidth = 0.4;
Figax = tight_subplot(1,2,gap,marg_h,marg_w);

% Make conc figure
axes(Figax(2))
new_S_vec = logspace(-8,-5,n_points_S*4);
P1P2_conc_vec = optimal_P1P2_concentrations(new_S_vec',alpha,alpha,beta,D_vals(1),Keq_vals(1),Keq_vals(1),gm);
P1C_conc_vec = optimal_P1C_concentrations(new_S_vec',alpha,alpha,beta,D_vals(1),Keq_vals(1),Keq_vals(1),gm);

physical_P1P2 = sum(P1P2_conc_vec > 0,2) >= 8;
physical_P1C = sum(P1C_conc_vec > 0,2) >= 7;

norm_vec = repmat(new_S_vec./D_vals(1),4,1);

colors1 = [64 176 166;64 176 166]/255;
colors1 = colors1.*[0.9;1.4];
dashtype = {'-',':','-',':'};

hold on
plot(new_S_vec(physical_P1P2),P1P2_conc_vec(physical_P1P2,4)'./norm_vec(4,physical_P1P2),dashtype{1},'LineWidth',LineWidth,'Color',colors1(1,:))
plot(new_S_vec(physical_P1P2),P1P2_conc_vec(physical_P1P2,end)'./norm_vec(end,physical_P1P2),dashtype{2},'LineWidth',LineWidth,'Color',colors1(1,:))
plot(new_S_vec(physical_P1C),P1C_conc_vec(physical_P1C,4)'./norm_vec(4,physical_P1C),dashtype{3},'LineWidth',LineWidth,'Color',colors1(2,:))
plot(new_S_vec(physical_P1C),P1C_conc_vec(physical_P1C,end)'./norm_vec(end,physical_P1C),dashtype{4},'LineWidth',LineWidth,'Color',colors1(2,:))

y_lim_vec = [0,0.7];
ylim(y_lim_vec)
yticks(y_lim_vec);
yticklabels({'0','0.7'})

xlim([new_S_vec(1),new_S_vec(end)])
xticks([new_S_vec(1),new_S_vec(end)]);
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

% Make panel B - minimum c1e ratio

colors = [220 38 127;...
    220 38 127]./255;
colors(2,:) = colors(2,:)*0.6;

axes(Figax(1))
hold on
D_labels = {'10^{-2}','10^{-3}'};
dashes = {'-',':'};
plot([1e-8,1e-5],[1,1],':','LineWidth',LineWidth,'color',0.5*[1,1,1])
for i = 1:length(Keq_vals)
    for j = 1:length(D_vals)
       %P1C_c1e = cellfun(@(x) x(4), all_c1e(i,j,:));
       %P1P2_c1e = cellfun(@(x) x(3), all_c1e(i,j,:));
       c1e_ratio = squeeze(P1C_P1P2_ratio_analytical(i,j,:));
       %c1e_ratio = P1C_c1e./P1P2_c1e;       
       plot(S_vec,squeeze(c1e_ratio),dashes{j},'Color',colors(i,:),'LineWidth',LineWidth)
    end
end


xlim([S_vec(1),S_vec(end)])
ylim([0.95,1.3])
xticks([S_vec(1),S_vec(end)]);
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

print(gcf, '-dpng','figs/SI_fig_BC_analytical_c0e_ratio.png','-r1200','-painters');

close all