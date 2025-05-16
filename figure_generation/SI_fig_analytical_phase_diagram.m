%% This script generates the SI figure with a phase diagram computed using
%the analytical approximation of the linear thermodynamic model

% This first portion of the script generates the phase boundary information
% for later plotting.
%
% The phase boundary data has been pre-generated, so you can skip
% running this portion.

clear;clc

%Set up parameters
Vc = 2e-15;
alpha= 1.8e-1;
beta = 2e2;
D_range = [1e-6,1e-3];
S_range = [1e-8,1e-5];
Keq_val = [1,10];
gm = 5e-4;
diagram_type = 'analytical';
n_points_S = 60;
n_bisections = 10;

%Loop through values of Keq
for i = 1:length(Keq_val)

    Keq = Keq_val(i);

    [upper_bounds,lower_bounds,phase_types,transition_points_cell,transition_map_cell]  ...
        = make_linear_phase_diagram(S_range,n_points_S,n_bisections,...
        alpha,alpha,beta,D_range,Vc,Keq,Keq,gm);

    save(['fig_data/phase_diagram_data_no_osm_Keq_',num2str(Keq),...
        '_',diagram_type,'.mat']);

end


%% Load data and make the phase diagrams

clear;clc

diagram_type = 'analytical'; 

newfigure(3.42/2, (1.4/3)*3.42/2);
gap = [0.2,0.15]; %height, width
marg_h = [0.2,0.15]; %lower upper
marg_w = [0.12,0.05]; %left right

Fig2ax = tight_subplot(1,2,gap,marg_h,marg_w);

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
    axes(Fig2ax(axis_inds(i)));
    hold on
    filename = ['fig_data/phase_diagram_data_no_osm_Keq_',num2str(Keq_val(i)),'_',diagram_type,'.mat'];
    load(filename);
    S_vec = logspace(log10(S_range(1)),log10(S_range(2)),n_points_S);
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
axes(Fig2ax(1));
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

print(gcf, '-dpng','figs/SI_fig_analytical_phase_diagram.png','-r1000','-painters');
