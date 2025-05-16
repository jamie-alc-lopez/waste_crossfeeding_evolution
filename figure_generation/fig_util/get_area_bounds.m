function [lower_bounds,upper_bounds] = get_area_bounds(phase_types,transition_points_cell,transition_map_cell,D_range)

%This script computes the upper and lower bounds of different phases
n_points = length(transition_points_cell);
lower_bounds = nan(length(phase_types),n_points);
upper_bounds = lower_bounds;

for i = 1:n_points
    transition_map = transition_map_cell{i};
    transition_points = transition_points_cell{i};    
    n_transitions = length(transition_points);
    
    %Look at first transition in map
    state = transition_map(1,1);
    matching_state = find(phase_types == state);
    lower_bounds(matching_state,i) = D_range(1);
    
    if n_transitions > 0
        upper_bounds(matching_state,i) = transition_points(1);
        
        for j = 1:(n_transitions-1)
            state = transition_map(j,2);
            matching_state = find(phase_types == state);
            lower_bounds(matching_state,i) = transition_points(j);
            upper_bounds(matching_state,i) = transition_points(j+1);
        end
        
        state = transition_map(end,2);
        matching_state = find(phase_types == state);
        lower_bounds(matching_state,i) = transition_points(end);
        upper_bounds(matching_state,i) = D_range(2);
        
    else
        upper_bounds(matching_state,i) = D_range(2);
    end
    
end

end

