%% This script generates the SI figure showing behaviors around the flux 
% plateau transition

clear;clc

% Plot the different growth functions
newfigure(3.42/2, (4.5/3)*3.42/2);
gap = [0.13,0.17]; %height, width
marg_h = [0.1,0.04]; %lower upper
marg_w = [0.145,0.13]; %left right
FontSize = 4;

E_colors = [254 97 0; 100 143 255; 120 94 240; 220 38 127; 255 176 0]/255; 
c_color = [86 180 233]/255;
dmat = [2 0 2 0; 1 1 1 1; 1/4 1/4 1/4 1/4];

letter_offsets = [-0.3,-0.25];

Figax = tight_subplot(3,2,gap,marg_h,marg_w);

%% Plot the P1 enzyme abundances
sim_file = 'fig_data/osm_sweep_Keq_1_no_sat_simulation_results.mat';

[lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file);
D_range = [sr.D_range(1),sr.D_range(end)];
S_range = [sr.supply_range(1),sr.supply_range(end)];

axes(Figax(1));
pos = get(gca,'Position');
pos(2) = pos(2)- 0.03;
set(gca,'Position',pos)

ifn = @(x) x(:,22); %57 for old fig, 22 for new fig
condition = 'D = 6.99e-6';
par_vec = ifn(sr.supply_mat);
LineWidth = 0.4;
xlab = '$\log_{10}$(Supply rate)';
tpt = 153;
tpt1 = mean(par_vec(tpt:(tpt+1)));

hold on

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
set(gca,'TickLabelInterpreter','latex')

plot([tpt1,tpt1],[1e-6,1e15],':','Color',[0.7,0.7,0.7])
P1_E1 = ifn(sr.P1_E_cell{1});
P1_T0 = ifn(sr.P1_E_cell{3});
P1_T1 = ifn(sr.P1_E_cell{4});

dashline(par_vec(1:end),P1_E1(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',E_colors(1,:));
dashline(par_vec(1:end),P1_T0(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',E_colors(3,:));
dashline(par_vec(1:end),P1_T1(1:end),dmat(3,1),dmat(3,2),dmat(3,3),dmat(3,4),'Color',E_colors(4,:));

%xlabel(xlab,'Interpreter','latex','FontSize',FontSize)
ylim([1e-6,1]); %6e11
yticks([1e-6,1e-3,1]);
yticklabels({'$10^{-6}$','$10^{-3}$','$10^0$'})
set(gca,'FontSize',4)
ylabel({'P1 enzymes'},'Interpreter','latex','FontSize',FontSize);
%yhandle = get(gca,'YLabel');
%pos = yhandle.Position;
%pos(1) = 0.6e-8;
%set(yhandle,'Position',pos);
set(gca,'YScale','log')
text(letter_offsets(1),1.12,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)

%% Plot the C enzyme abundances

axes(Figax(2));
pos = get(gca,'Position');
pos(2) = pos(2)- 0.03;
set(gca,'Position',pos)

hold on

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
set(gca,'TickLabelInterpreter','latex')

plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
C_E2 = ifn(sr.C_E_cell{2});
C_T1 = ifn(sr.C_E_cell{4});
C_T2 = ifn(sr.C_E_cell{5});

dashline(par_vec(1:end),C_E2(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',E_colors(2,:));
dashline(par_vec(1:end),C_T1(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',E_colors(5,:));
dashline(par_vec(1:end),C_T2(1:end),dmat(3,1),dmat(3,2),dmat(3,3),dmat(3,4),'Color',E_colors(4,:));

ylim([-0.1,1]); 
yticks([0,1]);
yticklabels({'$0$','$1$'})
set(gca,'FontSize',4)
ylabel({'C enzymes'},'Interpreter','latex','FontSize',FontSize);
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 0.6e-8;
set(yhandle,'Position',pos);

text(letter_offsets(2),1.12,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

leg_top = 1;
leg_left = 2e-8;
leg_colors =E_colors;
labels = {'$E_1$','$E_2$','$T_0$','$T_1$','$T_2$'};
line_length = 0.4e-7;
dashtype = {'-','-','-',':','-'};
spacing_x = 1e-8;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],dashtype{k},'Color',leg_colors(k,:))
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.11;
end


%% Plot the P1 internal concentrations

axes(Figax(3));

hold on

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
set(gca,'TickLabelInterpreter','latex')

plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
P1_c0 = ifn(sr.P1_c_cell{1});
P1_c1 = ifn(sr.P1_c_cell{2});

dashline(par_vec(1:end),P1_c0(1:end)-P1_c1(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',c_color);
%dashline(par_vec(1:end),P1_c1(1:end),dmat(2,1),dmat(2,2),dmat(2,3),dmat(2,4),'Color',c_color);

xlabel(xlab,'Interpreter','latex','FontSize',FontSize)
ylim([0,1e-2]);
yticks([0,1e-2])
yticklabels({'$0$','$10^{-2}$'})
set(gca,'FontSize',4)
ylabel({'P1 driving force, ($c_0-c_1$)'},'Interpreter','latex','FontSize',FontSize);

text(letter_offsets(1),1.14,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)


%% Plot the C internal concentrations

axes(Figax(4));

hold on

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
set(gca,'TickLabelInterpreter','latex')

plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
C_c1 = ifn(sr.C_c_cell{2});
C_c2 = ifn(sr.C_c_cell{3});

dashline(par_vec(1:end),C_c1(1:end)-C_c2(1:end),dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',c_color);
%dashline(par_vec(1:end),P1_c1(1:end),dmat(2,1),dmat(2,2),dmat(2,3),dmat(2,4),'Color',c_color);

xlabel(xlab,'Interpreter','latex','FontSize',FontSize)
ylim([0,1e-2]);
yticks([0,1e-2])
yticklabels({'$0$','$10^{-2}$'})
set(gca,'FontSize',4)
ylabel({'P2 driving force, ($c_1-c_2$)'},'Interpreter','latex','FontSize',FontSize);

text(letter_offsets(1),1.14,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)


%% Plot the overall biomass

axes(Figax(6))
axis off

axes(Figax(5));
pos = get(gca,'Position');
pos(1) = pos(1) + 0.8*pos(3);
set(gca,'Position',pos)
hold on

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
set(gca,'TickLabelInterpreter','latex')

plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
total_biomass = ifn(sr.P1_mat) + ifn(sr.P2_mat) + ifn(sr.C_mat);

plot(par_vec(1:end),total_biomass,'k-')

xlabel(xlab,'Interpreter','latex','FontSize',FontSize)
ylim([0,5e11]);
yticks([0,5e11])
yticklabels({'$0$','$5 \times 10^{11}$'})
set(gca,'FontSize',4)
ylabel({'Total community biomass,', '$\rho_\textrm{tot}$ (cells/L)'},'Interpreter','latex','FontSize',FontSize);

text(letter_offsets(1),1.12,'\textbf{E}','Interpreter','latex','Units','normalized','FontSize',4)


print(gcf, '-dpng','figs/SI_fig_flux_plateau.png','-r1200','-painters');
