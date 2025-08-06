%This script generates the SI figure analyzing the c0e vs. dilution
%relationship in a non-evolutionary condition and a condition with no
%thermodynamic inhibition

clear;clc

%Generate the blank figure
newfigure(3.42/2, (1.4/3)*3.42/2);
gap = [0.2,0.15]; %height, width
marg_h = [0.2,0.15]; %lower upper
marg_w = [0.17,0.05]; %left right

FigAx = tight_subplot(1,2,gap,marg_h,marg_w);



%% Simulate a non-evolutionary model under standard, no thermo, and no 
% thermo + no maintenance conditions

global par

%Load default parameter structure and set up parameters
load('../model_code/parameters/default_parameter_struct.mat');

n = 40;
D_vec = logspace(-6,-3,n);
par.m = 1;
par.supply = 3.107866187782010e-07; %Same as in Fig. 2D
E1P1 = 0.876100656900705; %Strategy from Fig. 2D
T0P1 = 0.061949671549648;
T1P1 = 0.061949671549648;
par.E = [E1P1 0];
par.T = [T0P1 T1P1 0];
min_growth = par.min_growth;
simu.rho = 1e10;
simu.c = zeros(2,3);

%Loop through and compute solution
sol = zeros(n,3);
for i = 1:n

    %Set dilution
    par.D = D_vec(i);

    %Simulation with thermo inhibition
    par.no_thermo = 0;
    par.min_growth = min_growth;
    [c_final,rho_final] = simulate_crossfeeding_ode15s(simu);

    if rho_final > 0
        c0e_standard(i) = c_final(end,1);
    else
        c0e_standard(i) = par.supply/par.D;
    end

    %Simulation without thermo inhibition
    par.no_thermo = 1;
    [c_final,rho_final] = simulate_crossfeeding_ode15s(simu);

    if rho_final > 0
        c0e_no_thermo(i) = c_final(end,1);
    else
        c0e_no_thermo(i) = par.supply/par.D;
    end

    %Simulation without thermo inhibition and without maintenance growth
    par.min_growth = 0;
    [c_final,rho_final] = simulate_crossfeeding_ode15s(simu);

    if rho_final > 0
        c0e_no_thermo_no_maint(i) = c_final(end,1);
    else
        c0e_no_thermo_no_maint(i) = par.supply/par.D;
    end
    

end

%% Plot the data

%Set up plotting variables
ext_c_color = [86 180 233]/255;
FontSize = 4;
dmat = [2 0 2 0; 1 1 1 1; 1/4 1/4 1/4 1/4];

ylimvec = [1e-6,10^(-1)];
ytickvec = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
ytickcell = {'$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'};

xlimvec = [1e-6,1e-3];
xtickvec = [1e-6,1e-5,1e-4,1e-3];
xtickcell = {'$-6$','$-5$','$-4$','$-3$'};


% (A), looking at simulation with thermodynamic inhibition
axes(FigAx(1));
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'YMinorTick','off')
set(gca,'XMinorTick','off')
hold on
dashline(D_vec,c0e_standard,dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',ext_c_color)
text(-0.3,1.08,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
xlabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
ylabel({'Extracellular', 'compound 0'},'Interpreter','latex')
set(gca,'FontSize',FontSize)
xlim(xlimvec)
xticks(xtickvec);
xticklabels(xtickcell)
ylim(ylimvec)
yticks(ytickvec);
yticklabels(ytickcell)
set(gca,'TickLabelInterpreter','latex')


% (B), looking at simulation with no thermodynamic inhibition and a
% simulation with no thermo + no maintenance growth rate
axes(FigAx(2));
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'YMinorTick','off')
set(gca,'XMinorTick','off')
hold on
dashline(D_vec,c0e_no_thermo,dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',ext_c_color)
dashline(D_vec,c0e_no_thermo_no_maint,dmat(2,1),dmat(2,2),dmat(2,3),dmat(2,4),'Color','k')
text(-0.3,1.08,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)
xlabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
set(gca,'FontSize',FontSize)
xlim(xlimvec)
xticks(xtickvec);
xticklabels(xtickcell)
ylim(ylimvec)
yticks(ytickvec);
yticklabels(ytickcell)
set(gca,'TickLabelInterpreter','latex')


%Legend
colors = {ext_c_color,'k'};
leg_top = 1e-1;
leg_left = 0.5e-5;
leg_labels_c = {'$g_\textrm{m}= 5 \times 10^{-4}$','$g_\textrm{m} = 0$'};
line_length = 15e-6;
spacing_x = 0.5e-5;
curr_y = leg_top;
for k = 1:length(leg_labels_c)
    dashline([leg_left,leg_left + line_length],[curr_y,curr_y],...
        dmat(k,1),dmat(k,2),dmat(k,3),dmat(k,4),'Color',colors{k})
    text(leg_left + line_length + spacing_x,curr_y,leg_labels_c{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y*0.2;
end



print(gcf, '-dpng','figs/SI_fig_c0e_no_inhibition.png','-r800','-painters');
