%% This script generates the SI fig showing bistability as a function of 
%initial P1 abundance

clear;clc

%Set up simulation parameters
global par
load('../model_code/parameters/default_parameter_struct.mat')

%Particular parameters for bistable example simulation. Parameters taken
%from an evolution simulation that initially showed this behavior
par.osm_coeff = 1;
par.E = [0.872833717775150,0;0,0.877896422740764];
par.T = [0.0605001503358823,0.0666661318889683,0;0,0.0417604005370844,0.0803431767221518];
par.D = 1e-6;
par. m = 2;
par.supply = 1.482020705798860e-6;
n = 20;
rho1_vec = logspace(8,11,n);
rho1_vec(1) = 1;
rho2 = 705166807262.602;
simu.c = [0.493634155293565,0.492459874970296,0;4.10581669642876e-22,0.492333544654768,0.491935684019434;0.493724268614594,0.492378096272252,0.491912592351409];
par.plt = 0;
par.no_T1 = 0;
par.tol = 1e-13;

%Run simulation across different initial conditions
for i = 1:n
    simu.rho = [rho1_vec(i);rho2];
    [c_final,rho_final,c_store{i},rho_store{i},t_store{i}] = simulate_crossfeeding_ode15s(simu);
end


%Plot the population dynamics across initial conditions
figure
pos = get(gcf,'Position');
pos(3) = 1.5*pos(3);
pos(4) = 1.3*pos(4);
set(gcf,'Position',pos)

sim_init = 8;
sim_increment = 1;
plt_range = sim_init:sim_increment:n;
c_offset = 1;
colors = parula(length(plt_range)+c_offset);
colors = colors(1:end-c_offset,:);
pop_ylim_vec = [1e6,1e13];
pop_ytick_vec = [1e6,1e8,1e10,1e12];
time_lim = [1e0,1e8]; %[1e1,21560600];
time_tick_vec = [1e0,1e2,1e4,1e6,1e8];
LineWidth = 1.5;
FontSize = 17;

%Plot the P1 dynamics
subplot(2,2,1)
hold on
tn = 1;
for i = plt_range
    plot(t_store{i},rho_store{i}(1,:),'-','LineWidth',LineWidth,...
        'Color',colors(tn,:))
    tn = tn + 1;
end

set(gca,'XScale','log')
set(gca,'YScale','log')
ylabel('P1 population (cells/L)','Interpreter','latex')
ylim(pop_ylim_vec)
yticks(pop_ytick_vec)
xlim(time_lim)
xticks(time_tick_vec)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)
text(-0.1,1.1,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',FontSize)


%Plot the P2 dynamics
subplot(2,2,2)
hold on
tn = 1;
for i = plt_range
    plot(t_store{i},rho_store{i}(2,:),'-','LineWidth',LineWidth,...
        'Color',colors(tn,:))
    tn = tn + 1;
end
set(gca,'XScale','log')
set(gca,'YScale','log')
ylabel('P2 population (cells/L)','Interpreter','latex')
ylim(pop_ylim_vec)
yticks(pop_ytick_vec)
xlim(time_lim)
xticks(time_tick_vec)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)
text(-0.1,1.1,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',FontSize)


%Plot the c0e dynamics
subplot(2,2,3)
hold on
tn = 1;
for i = plt_range
    plot(t_store{i},squeeze(c_store{i}(end,1,:)),'-','LineWidth',LineWidth,...
        'Color',colors(tn,:))
    tn = tn + 1;
end
set(gca,'XScale','log')
%set(gca,'YScale','log')
xlabel('Time','Interpreter','latex')
ylabel('Extracellular compound 0','Interpreter','latex')
ylim([0,2])
xlim(time_lim)
xticks(time_tick_vec)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)
text(-0.1,1.1,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',FontSize)


%Plot the c2e dynamics
subplot(2,2,4)
hold on
tn = 1;
for i = plt_range
    plot(t_store{i},squeeze(c_store{i}(end,3,:)),'-','LineWidth',LineWidth,...
        'Color',colors(tn,:))
    tn = tn + 1;
end
set(gca,'XScale','log')
%set(gca,'YScale','log')
xlabel('Time','Interpreter','latex')
ylabel('Extracellular compound 2','Interpreter','latex')
ylim([0,0.6])
yticks([0,0.2,0.4,0.6])
xlim(time_lim)
xticks(time_tick_vec)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)
text(-0.1,1.1,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',FontSize)


print(gcf, '-dpng','figs/SI_fig_bistability.png','-r800','-painters');



