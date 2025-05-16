function [lower_bounds,upper_bounds,phase_types,sr] = make_osm_phase_diagram(sim_file)

%This script identifies phase boundaries based on numerical simulation results

%Load in processed simulation data
load(sim_file);

%Loop through supply space and identify transition points + class in
%dilution space
S_vec = sr.supply_range;
D_vec = sr.D_range;
D_range = [sr.D_range(1),sr.D_range(end)];
n_points = length(S_vec);
phase_matrix = sr.phase_matrix;
for i = 1:n_points
    phase_vec = phase_matrix(i,:);
    transition_ind = find(diff(phase_vec));
    n_transitions = length(transition_ind);
    transition_points_cell{i} = nan(n_transitions,1);
    transition_map_cell{i} = nan(n_transitions,2);
    
    %Identify each transition
    for j = 1:n_transitions
        transition_points_cell{i}(j) = mean(D_vec(transition_ind(j):(transition_ind(j)+1)));
        transition_map_cell{i}(j,1) = phase_vec(transition_ind(j));
        transition_map_cell{i}(j,2) = phase_vec(transition_ind(j)+1);
    end
    
    if n_transitions == 0
        transition_map_cell{i}(1,1) = phase_vec(1);
        transition_map_cell{i}(1,2) = phase_vec(2);
    end
    
end

% %Sort transition points to go from least to greatest
% for i = 1:length(transition_map_cell)
%     transition_points = transition_points_cell{i};
%     transition_map = transition_map_cell{i};
%     [transition_points,ia] = sort(transition_points);
%     transition_points_cell{i} = transition_points;
%     transition_map_cell{i} = transition_map(ia,:);
% end

%Get phase area bounds
phase_types = unique(phase_matrix(:));
[lower_bounds,upper_bounds] = get_area_bounds(phase_types,transition_points_cell,transition_map_cell,D_range);


end

