%% This script generates fig 4 of the manuscript

% Plot the different growth functions
newfigure(3.42/2, (4.0/3)*3.42/2);
gap = [0.13,0.17]; %height, width
marg_h = [0.1,0.04]; %lower upper
marg_w = [0.145,0.13]; %left right
FontSize = 4;

Fig4ax = tight_subplot(3,2,gap,marg_h,marg_w);

axes(Fig4ax(2))
axis off

axes(Fig4ax(1))
text(-0.35,1.1,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
hold on

pos = get(gca,'Position');
pos(1) = 2.45*pos(1);
set(gca,'Position',pos)

non_osm_fun = @(c1,c2,Keq) max([c1 - c2./Keq; zeros(size(c1))],[],1);
osm_fun = @(c1,c2,Keq) (1 - c1 - c2).*max([c1 - c2./Keq; zeros(size(c1))],[],1);

c1_vec = linspace(0,1,100);
c2_val = 0;
Keq = 2;
color = 'k';
plot(c1_vec,osm_fun(c1_vec,c2_val,Keq),'-','Color',color)
plot(c1_vec,non_osm_fun(c1_vec,c2_val,Keq),':','Color',color)


xlim([0,1])
ylim([0,1.2])
xticks([])
yticks([])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FontSize)
ylabel('Growth rate','Interpreter','latex','FontSize',FontSize)
xlabel({'Internal compound', 'concentration'},'Interpreter','latex','FontSize',FontSize)

leg_top = 1.2;
leg_left = 7e-2;
mylines = {':','-'};
leg_colors = {'k','k'};
labels = {'Linear model','Osmotic model'};
line_length = 15e-2;
spacing_x = 0.05;
curr_y = leg_top;
for k = 1:length(labels)
    plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',leg_colors{k})
    text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
    curr_y = curr_y-0.12;
end


%% Make osmotic phase diagrams

Keq_val = [1,10];
axis_inds = [3,4];
labels = {'\textbf{B}','\textbf{C}'};
colors = [255 255 255; ...
    240 228 66;...
    213 94 0;...
    204 102 119;...
    0 158 115;...
    100 100 100];
%colors order is null, P1, P2, P1P2, P1C
colors = colors./255;
colors(2,:) = 0.8*colors(2,:);
letter_offsets = [-0.35,-0.25];
highlighted_D = -5.6834;

for i = 1:length(Keq_val)

    Keq = Keq_val(i);
    sim_file = ['fig_data/osm_sweep_Keq_',num2str(Keq),'_no_sat_simulation_results.mat'];
    [lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file);
    phase_types(phase_types == -1) = 5;

    D_range = [sr.D_range(1),sr.D_range(end)];
    S_range = [sr.supply_range(1),sr.supply_range(end)];

    axes(Fig4ax(axis_inds(i)));
    hold on
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

    xhandle = get(gca,'XLabel');
    pos = xhandle.Position;
    pos(2) = 0.6*pos(2);
    set(xhandle,'Position',pos);

    if i == 1
        yhandle = get(gca,'YLabel');
        pos = yhandle.Position;
        %pos(1) = 0.1*pos(1);
        set(yhandle,'Position',pos);

        plot(log10([S_range(1) S_range(2)]),[highlighted_D,highlighted_D],'k:')
    end

    xhandle = get(gca,'XLabel');
    pos = xhandle.Position;
    pos(2) = -6.4;
    set(xhandle,'Position',pos);

end

%Make legend
axes(Fig4ax(3));
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
point1 = [-8.1,level];
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

point5 = point4 + [1.5*point_movement,0];
make_rec(point5,colors(5,:))
text(point5(1) + point_offset,point5(2),'P1+C',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

point6 = point5 + [1.4*point_movement,0];
make_rec(point6,colors(6,:))
text(point6(1) + point_offset,point6(2),'Osc.',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',4);

%% Show discontinuous properties across the cross-feeding transition

sim_file = 'fig_data/osm_sweep_Keq_1_no_sat_simulation_results.mat';
[lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file);
D_range = [sr.D_range(1),sr.D_range(end)];
S_range = [sr.supply_range(1),sr.supply_range(end)];

axes(Fig4ax(5));
pos = get(gca,'Position');
pos(2) = pos(2)- 0.03;
set(gca,'Position',pos)

ifn = @(x) x(:,22); %57 for old fig, 22 for new fig
condition = 'D = 6.99e-6';
par_vec = ifn(sr.supply_mat);
LineWidth = 0.4;
tpt = 153;
tpt1 = mean(par_vec(tpt:(tpt+1)));

hold on
colorset = [0.5*colors(2,:); 250 50 32];
colorset(2,:) = 0.9*colorset(2,:)/255;
colororder(colorset)

%Left axis, biomass
yyaxis left
plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
P1_biomass = ifn(sr.P1_mat);%+ifn(sr.P2_mat)+ifn(sr.C_mat);
plot(par_vec(1:end),P1_biomass(1:end),'-','LineWidth',LineWidth,'Color',2*colorset(1,:));
%plot(par_vec((tpt+1):end),total_biomass((tpt+1):end),'-','LineWidth',LineWidth,'Color',colorset(1,:));
ylim([0,6e11]); %6e11
yticks([0,6e11]);
yticklabels({'$0$','$6 \times 10^{11}$'})
set(gca,'FontSize',4)
ylabel({'P1 population','(cells/L)'},'Interpreter','latex','FontSize',FontSize);
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 0.6e-8;
set(yhandle,'Position',pos);

%Right axis, osmotic burden
yyaxis right
burden = ifn(sr.burden_mat);
plot(par_vec(1:end),burden(1:end),'-','LineWidth',LineWidth,'Color',colorset(2,:));
%plot(par_vec((tpt+1):end),burden((tpt+1):end),'-','LineWidth',LineWidth,'Color',colorset(2,:));
yticks([0,1]);
ylim([0,1.1])
ylabel({'Mean toxic burden'},'Interpreter','latex','FontSize',FontSize);


%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
xlabel('$\log_{10}$(Supply rate)','Interpreter','latex','FontSize',FontSize)
set(gca,'TickLabelInterpreter','latex')

text(letter_offsets(1),1.12,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)


%% Show continuous properties across the cross-feeding transition

axes(Fig4ax(6))
pos = get(gca,'Position');
pos(2) = pos(2)- 0.03;
set(gca,'Position',pos)
colorset = [86 180 233; 93 58 155]/255;
colorset(1,:) = 0.7*colorset(1,:);
colororder(colorset)

%Right axis, total flux
hold on
yyaxis right
plot([tpt1,tpt1],[-50,1e15],':','Color',[0.7,0.7,0.7])
plot(par_vec,ifn(sr.P1_flux_mat.*sr.P1_mat) + ifn(sr.P2_flux_mat.*sr.P2_mat),'-','LineWidth',LineWidth,'Color',colorset(2,:));
ylims = get(gca,'YLim');
ylim([0,0.5e-5]);
yticks([0,0.5e-5]);
yticklabels({'$0$','$5 \times 10^{-6}$'})
set(gca,'FontSize',4)
ylabel('Total flux','Interpreter','latex','FontSize',FontSize);
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 0.18e-4;
set(yhandle,'Position',pos);
set(gca,'Clipping','on')


%Left axis, extracellular c1
yyaxis left
plot(par_vec,log10(ifn(sr.extra_c_cell{1})),'-','LineWidth',LineWidth,'Color',colorset(1,:))
yticks([-4,1]);
ylim([-4,1])
yticklabels({-4,1})
ylabel({'$\log_{10}$($c_{0,e}$)'},'Interpreter','latex','FontSize',FontSize);
yhandle = get(gca,'YLabel');
pos = yhandle.Position;
pos(1) = 0.5e-8;
set(yhandle,'Position',pos);

%Modify x axis
set(gca,'XScale','log')
xlim([S_range(1),S_range(2)]);
xticks([S_range(1),S_range(2)]);
xticklabels({'$-8$','$-5$'})
xlabel('$\log_{10}$(Supply rate)','Interpreter','latex','FontSize',FontSize)
set(gca,'TickLabelInterpreter','latex')

set(gca,'Clipping','on')

text(letter_offsets(1),1.12,'\textbf{E}','Interpreter','latex','Units','normalized','FontSize',4)

print(gcf, '-dpng','figs/fig4.png','-r1200','-painters');


