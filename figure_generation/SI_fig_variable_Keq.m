%% This script generates the phase diagram for the case of no energy yield
% from the second reaction
% First section: generate phase boundaries from the numerical simulation results

clear;clc

D_range = [1e-6,1e-3];
S_range = [1e-8,1e-5]; 

sim_file = ['fig_data/no_osm_sweep_Keq1_1_Keq2_0p1_no_sat_simulation_results.mat'];
[lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file);
save(['fig_data/phase_diagram_data_no_osm_sweep_Keq1_1_Keq2_0p1_no_sat_simulation_results.mat']);


%% Make the phase diagrams

clear;clc

newfigure(0.9*3.42/2, 0.7*3.42/2);

pos = get(gca,'Position');
pos(4) = 0.94*pos(4);
pos(2) = 1.3*pos(2);
pos(3) = 0.93*pos(3);
pos(1) = 1.4*pos(1);
set(gca,'Position',pos)


colors = [255 255 255; ...
    240 228 66;...
    213 94 0;...
    204 102 119;...
    0 158 115];
%colors order is null, P1, P2, P1P2, P1C
colors = colors./255;
colors(2,:) = 0.8*colors(2,:);
FontSize = 4;

hold on 
filename = ['fig_data/phase_diagram_data_no_osm_sweep_Keq1_1_Keq2_0p1_no_sat_simulation_results.mat'];
load(filename);
S_vec = sr.supply_range';
for j = 1:length(phase_types)
    nan_ind = ~isnan(lower_bounds(j,:));
    x = [S_vec(nan_ind),fliplr(S_vec(nan_ind))];
    y = [lower_bounds(j,nan_ind),fliplr(upper_bounds(j,nan_ind))];
    fill(log10(x),log10(y),colors(phase_types(j)+1,:),'LineWidth',0.6,'EdgeColor','k')
end
xlim(log10(S_range));
ylim(log10(D_range));
box on
set(gca,'FontSize',4)

ylabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
xlabel('$\log_{10}$(Supply rate)','Interpreter','latex','FontSize',FontSize)

xticks(log10([S_range(1),S_range(2)]));
yticks(log10([D_range(1),D_range(2)]));
xticklabels({'$-8$','$-5$'})
yticklabels({'$-6$','$-3$'})
set(gca,'TickLabelInterpreter','latex')

% xhandle = get(gca,'XLabel');
% pos = xhandle.Position;
% pos(2) = 0.6*pos(2);
% set(xhandle,'Position',pos);

%Make legend
ax = gca;
ax.Clipping = 'off';
legx = log10(S_range(2));
level = -6.8;
rec_width = abs(0.11*1.5);
rec_height = rec_width;
point_movement = -0.13*legx;
make_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/2,coord(2)-0.92*rec_height/2,rec_width,rec_height],...
    'FaceColor',color,'EdgeColor','k','LineWidth',0.6);

point_offset = -0.03*legx;
point1 = [-7.65,level];
make_rec(point1,colors(1,:))
text(point1(1) + point_offset,point1(2),'$\emptyset$',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point2 = point1 + [0.8*point_movement,0];
make_rec(point2,colors(2,:))
text(point2(1) + point_offset,point2(2),'P1',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point3 = point2 + [point_movement,0];
make_rec(point3,colors(3,:))
text(point3(1) + point_offset,point3(2),'P2',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point4 = point3 + [point_movement,0];
make_rec(point4,colors(4,:))
text(point4(1) + point_offset,point4(2),'P1+P2',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

print(gcf, '-dpng','figs/SI_fig_variable_Keq.png','-r800','-painters');


