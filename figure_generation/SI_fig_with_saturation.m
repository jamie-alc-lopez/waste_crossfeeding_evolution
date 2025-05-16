%% This script generates the SI figure showing a version of Figure 2 but with
% a saturating growth function

% This first portion of the script generates the phase boundary information
% for later plotting.
%
% The phase boundary data has been pre-generated, so you can skip
% running this portion.

clear;clc

%Set up parameters
D_range = [1e-6,1e-3];
S_range = [1e-8,1e-5];
Keq_val = [1,10];

%Loop through values of Keq
for i = 1:length(Keq_val)

    Keq = Keq_val(i);

    %Purely numerical phase diagram from pre-run full model numerics
    n_points_S = 200;
    n_points_D = 200;
    sim_file = ['fig_data/no_osm_sweep_Keq_',num2str(Keq),...
        '_sat_simulation_results.mat'];
    [lower_bounds,upper_bounds,phase_types,sr] = ...
        make_osm_phase_diagram(sim_file);

    save(['fig_data/phase_diagram_data_no_osm_Keq_',num2str(Keq),...
        '_with_saturation','.mat']);

end



%% Load data and make the phase diagrams

clear;clc

highlighted_S = 3.1079e-07;

newfigure(3.42/2, (3/3)*3.42/2);
gap = [0.2,0.15]; %height, width
marg_h = [0.10,0.08]; %lower upper
marg_w = [0.12,0.05]; %left right

SIFigax = tight_subplot(2,2,gap,marg_h,marg_w);

Keq_val = [1,10];
axis_inds = [1,2];
labels = {'\textbf{A}','\textbf{B}','\textbf{C}','\textbf{D}'};
colors = [255 255 255; ...
    240 228 66;...
    213 94 0;...
    204 102 119;...
    0 158 115];
%color order is null, P1, P2, P1P2, P1C
colors = colors./255;
colors(2,:) = 0.8*colors(2,:);
FontSize = 4;
letter_offsets = [-0.35,-0.25,-0.35,-0.25];

for i = 1:length(Keq_val)
    axes(SIFigax(axis_inds(i)));
    hold on
    filename = ['fig_data/phase_diagram_data_no_osm_Keq_',...
        num2str(Keq_val(i)),'_with_saturation','.mat'];
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

    if i == 1
        ylabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
        plot(log10([highlighted_S,highlighted_S]),log10([D_range(1) D_range(2)]),'k:')

    end
    xlabel('$\log_{10}$(Supply rate)','Interpreter','latex','FontSize',FontSize)

    text(letter_offsets(i),1.12,labels{i},'Interpreter','latex','Units','normalized','FontSize',4)
    figlabel = ['$K_{\textrm{eq}} = ',num2str(Keq),'$'];
    text(0.5,1.1,figlabel,'Interpreter','latex','Units',...
        'normalized','FontSize',4,'HorizontalAlignment','center')

    xticks(log10([S_range(1),S_range(2)]));
    yticks(log10([D_range(1),D_range(2)]));
    xticklabels({'$-8$','$-5$'})
    yticklabels({'$-6$','$-3$'})
    set(gca,'TickLabelInterpreter','latex')

    if i == 1
        yhandle = get(gca,'YLabel');
        pos = yhandle.Position;
        %pos(1) = 0.1*pos(1);
        set(yhandle,'Position',pos);
    end

    xhandle = get(gca,'XLabel');
    pos = xhandle.Position;
    pos(2) = -6.4;
    set(xhandle,'Position',pos);
end

%Make legend
axes(SIFigax(1));
ax = gca;
ax.Clipping = 'off';
legx = log10(S_range(2));
level = -7.2;
rec_width = abs(0.11*3);
rec_height = rec_width;
point_movement = -0.26*legx;
make_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/2,coord(2)-0.92*rec_height/2,rec_width,rec_height],...
    'FaceColor',color,'EdgeColor','k','LineWidth',0.6);

point_offset = -0.06*legx;
point1 = [-6.5,level];
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


%% Plot populations and external c across dilution rates

load('fig_data/no_osm_sweep_Keq_1_sat_simulation_results.mat');

stripe_ind = 100;

ext_c_colors = [86 180 233]/255;
rho_colors = [68 170 153; 136 34 85]/255;
dmat = [2 0 2 0; 1 1 1 1; 1/4 1/4 1/4 1/4];

c1_vec = sr.extra_c_cell{1}(:,stripe_ind);
c2_vec = sr.extra_c_cell{2}(:,stripe_ind);
c3_vec = sr.extra_c_cell{3}(:,stripe_ind);

P1_vec = sr.P1_mat(:,stripe_ind);
P2_vec = sr.P2_mat(:,stripe_ind);

D_vec = flip(sr.D_range);

%External concentrations
axes(SIFigax(4));
hold on
dashline(D_vec,c1_vec,dmat(1,1),dmat(1,2),dmat(1,3),dmat(1,4),'Color',ext_c_colors)
dashline(D_vec,c2_vec,dmat(2,1),dmat(2,2),dmat(2,3),dmat(2,4),'Color',ext_c_colors)
dashline(D_vec,c3_vec,dmat(3,1),dmat(3,2),dmat(3,3),dmat(3,4),'Color',ext_c_colors)
text(letter_offsets(3),1.08,labels{4},'Interpreter','latex','Units','normalized','FontSize',4)
xlabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
ylabel('Concentration','Interpreter','latex')
set(gca,'FontSize',FontSize)
set(gca,'XScale','log')
xticks([D_range(1),D_range(2)]);
xlim([D_range(1),D_range(2)]);
xticklabels({'$-6$','$-3$'})
yticks([0,3e-2])
ylim([0,3e-2]);
yticklabels({'0','$3 \times 10^{-2}$'})
set(gca,'TickLabelInterpreter','latex')

leg_top = 3e-2;
leg_left = 0.5e-4;
leg_labels_c = {'$c_{0,e}$','$c_{1,e}$','$c_{2,e}$'};
line_length = 15e-5;
spacing_x = 0.5e-4;
curr_y = leg_top;
for k = 1:length(leg_labels_c)
    dashline([leg_left,leg_left + line_length],[curr_y,curr_y],...
        dmat(k,1),dmat(k,2),dmat(k,3),dmat(k,4),'Color',ext_c_colors)
    text(leg_left + line_length + spacing_x,curr_y,leg_labels_c{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.3e-2;
end

yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 0.6e-6;
set(yhandle,'Position',pos);


%Population densities
axes(SIFigax(3));
hold on
plot(D_vec,P1_vec,'-','Color',colors(2,:))
plot(D_vec,P2_vec,'-','Color',colors(3,:))
text(letter_offsets(4),1.08,labels{3},'Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'FontSize',FontSize)
xlabel('$\log_{10}$(Dilution rate)','Interpreter','latex','FontSize',FontSize)
ylabel({'Population (cells/L)'},'Interpreter','latex')
set(gca,'XScale','log')
xticks([D_range(1),D_range(2)]);
xlim([D_range(1),D_range(2)]);
xticklabels({'$-6$','$-3$'})
yticks([1e10,1e13])
ylim([10^(10),1e13]);
yticklabels({'$10^{10}$','$10^{13}$'})
set(gca,'YScale','log')
set(gca,'TickLabelInterpreter','latex')

leg_top = 1e13;
leg_left = 0.5e-4;
mylines = {'-','-','-'};
leg_labels_rho = {'P1','P2'};
line_length = 15e-5;
spacing_x = 0.5e-4;
spacing_y = 6e-2*1e14;
curr_y = leg_top;
for k = 1:length(leg_labels_rho)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k+1,:))
    text(leg_left + line_length + spacing_x,curr_y,leg_labels_rho{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-spacing_y;
end

print(gcf, '-dpng','figs/SI_fig_with_saturation.png','-r800','-painters');

