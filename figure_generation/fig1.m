%% This script generates the evolution trajectory depicted in Fig 1B. The 
%rest of the figure is generated via biorender

clear;clc

%First, generate the trajectory used in this figure
global par
load('../model_code/parameters/default_parameter_struct.mat');
par.supply = 1e-5;
par.D = 1e-4; 
par.osm_coeff = 0;

[current_ET,par_store, opt_growth,total_rho,species_enzymes,simu] = ...
    chemostat_invasion_ode15s();

save('fig_data/P1P2_evolution_trajectory.mat')

%% Plot the trajectories

colors = [255 255 255; ...
    240 228 66;...
    213 94 0;...
    204 102 119;...
    0 158 115];
%colors order is null, P1, P2, P1P2, P1C
colors = colors./255;
load('fig_data/P1P2_evolution_trajectory.mat')
species_enzymes = fliplr(species_enzymes);
color_ind = [2 5 3];

for i = 1:size(species_enzymes,1)
figure
barh(species_enzymes(i,:),0.7,'FaceColor',colors(color_ind(i),:),'LineWidth',1.5)
set(gca,'XScale','log')
xlim([1e-3,1])
box off
xticks([])
yticks([1 2 3 4 5])
yticklabels({})
ylim([0.5,5.5])
H = gca;
H.LineWidth = 1.5;
pos = get(gcf,'Position');
pos(3) = 0.15*pos(3);
pos(4) = 0.25*pos(4);
set(gcf,'Position',pos)
print(gcf, '-dpng',['figs/fig1B_P1P2_evolution_',num2str(i),'.png'],'-r1200');

end