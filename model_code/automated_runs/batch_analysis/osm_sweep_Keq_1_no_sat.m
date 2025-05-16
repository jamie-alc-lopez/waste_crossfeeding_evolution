clear;clc

%This script imports and processes the simulation data from the osmotic 
%model with Keq 1. It also plots a rough phase diagram for diagnostic
%purposes, the final phase diagram scripts are in the figure_generation
%folder

sim_name = 'osm_sweep_Keq_1_no_sat';
num_jobs = 800;
run_per_job = 50;
E_cutoff = 1e-6;
enzyme_only = 1;
[parameter_table,sr] = import_batch_run(sim_name,num_jobs,run_per_job,E_cutoff);

sr.burden_mat = zeros(size(sr.raw_c_cell));

for l = 1:size(sr.raw_c_cell,1)
    for m = 1:size(sr.raw_c_cell,2)
        rho = sr.raw_rho_cell{l,m};
        c = sr.raw_c_cell{l,m};
        burden = 0;
        for n = 1:length(rho)
            burden = burden + rho(n)*sum(c(n,:));
        end
        sr.burden_mat(l,m) = burden/sum(rho);    
    end
end
        
save(['../../../figure_generation/fig_data/',sim_name,'_simulation_results.mat'])

%% Plot phase diagram

%unique_phases = unique(sr.phase_matrix(~isnan(sr.phase_matrix(:))));
unique_phases = [-2 -1 0 1 2 3 4 5];
mymap = parula(length(unique_phases));
colormap(mymap)
data = flipud(sr.phase_matrix');
figure
imagesc(data,'AlphaData',~isnan(data),[-1,5])

xlabel('Nutrient supply, S','Interpreter','latex')
ylabel('Dilution rate, D','Interpreter','latex')
xticks([0.5,length(sr.supply_range)+ 0.5])
yticks([0.5,length(sr.supply_range)+ 0.5])
xticklabels({num2str(min(sr.supply_range)),num2str(max(sr.supply_range))})
yticklabels({num2str(max(sr.D_range)),num2str(min(sr.D_range))})
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)


%Make legend
ax = gca;
pos = ax.Position;
scale = 0.9;
pos(4) = scale*pos(4);
pos(2) = pos(2) + (1/scale-1)*pos(4);
ax.Position = pos;
ax.Clipping = 'off';

rec_width = 5;
make_rec = @(coord,color)rectangle('Position',...
    [coord(1)-rec_width/2,coord(2)-rec_width/2,rec_width,rec_width],...
    'FaceColor',color,'EdgeColor',color);

point_mult = 8;
point_offset = 10*0.5;
point1 = [0,170];
make_rec(point1,mymap(unique_phases==0,:))
text(point1(1) + point_offset,point1(2),'$\emptyset$',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',18);

point2 = point1 + point_mult*[3,0];
make_rec(point2,mymap(unique_phases==1,:))
text(point2(1) + point_offset,point2(2),'P1',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',14);

point3 = point2 + point_mult*[3.5,0];
make_rec(point3,mymap(unique_phases==2,:))
text(point3(1) + point_offset,point3(2),'P2',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',14);

point4 = point3 + point_mult*[3.5,0];
make_rec(point4,mymap(unique_phases==3,:))
text(point4(1) + point_offset,point4(2),'P1/P2',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',14);

point5 = point4 + point_mult*[4.5,0];
make_rec(point5,mymap(unique_phases==4,:))
text(point5(1) + point_offset,point5(2),'P1/C',...
    'HorizontalAlignment','left','VerticalAlignment','middle',...
    'Interpreter','latex','FontSize',14);
