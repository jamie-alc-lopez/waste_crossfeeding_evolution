function [upper_bounds,lower_bounds,phase_types,transition_points_cell,transition_map_cell] = ...
    make_linear_phase_diagram(S_range,n_points,n_bisections,alpha1,alpha2,...
    beta,D_range,Vc,Keq1,Keq2,gm)

%This script computes a phase diagram based on the analytical approximation
%of the linear thermodynamics model 

S_vec = logspace(log10(S_range(1)),log10(S_range(2)),n_points);

%Loop through grid
for i = 1:n_points
   S = S_vec(i); 
   [transition_points_cell{i},transition_map_cell{i}] = ...
      find_linear_phase_boundaries(S,n_points,n_bisections,alpha1,alpha2,...
      beta,D_range,Vc,Keq1,Keq2,gm); 
    
end

%Sort transition points to go from least to greatest
for i = 1:length(transition_map_cell)
    transition_points = transition_points_cell{i};
    transition_map = transition_map_cell{i};
    [transition_points,ia] = sort(transition_points);
    transition_points_cell{i} = transition_points;
    transition_map_cell{i} = transition_map(ia,:);
end


%Get phase boundary lines
transition_types = unique(cell2mat(transition_map_cell'),'rows');
transition_matrix = nan(size(transition_types,1),n_points);
for i = 1:size(transition_types,1)
    transition = transition_types(i,:);
    for j = 1:n_points
        [~,ia] = intersect(transition_map_cell{j},transition,'rows');
        if ~isempty(ia)
            transition_matrix(i,j) = transition_points_cell{j}(ia);
        end
    end
end

%Get phase transition lines
transition_types = unique(cell2mat(transition_map_cell'),'rows');
transition_matrix = nan(size(transition_types,1),n_points);
for i = 1:size(transition_types,1)
    transition = transition_types(i,:);
    for j = 1:n_points
        [~,ia] = intersect(transition_map_cell{j},transition,'rows');
        if ~isempty(ia)
            transition_matrix(i,j) = transition_points_cell{j}(ia);
        end
    end
end

%Get phase areas
phase_types = unique(transition_types(:));
lower_bounds = nan(length(phase_types),n_points);
upper_bounds = lower_bounds;
for i = 1:length(phase_types)
    phase_type = phase_types(i);
    for j = 1:n_points
        transition_map = transition_map_cell{j};
        transition_points = transition_points_cell{j};
        matching_states = transition_map == phase_type; 
        
        %Check whether the state occurs once or twice in the transition map
        if sum(matching_states(:)) == 2
            
            %If no transition occurs, make whole region single state
            if size(transition_map,1) == 1
               lower_bounds(i,j) = D_range(1);
               upper_bounds(i,j) = D_range(2);                 
               
            %If intermediate transition occurs, use transition points to
            %assign bounds
            else
               matching_transitions = find(sum(matching_states,2));
               lower_bounds(i,j) = transition_points(matching_transitions(1));
               upper_bounds(i,j) = transition_points(matching_transitions(2));     
            end
            
        elseif sum(matching_states(:)) == 1
            %Only appearing once implies a transition at the domain bounds
            matching_position = sum(matching_states,1);
            
            %Transition at lower domain bound
            if matching_position(1)
                lower_bounds(i,j) = D_range(1);
                upper_bounds(i,j) = transition_points(1);
                if j > 1
                    if isnan(lower_bounds(i,j-1))
                        lower_bounds(i,j-1) = D_range(1);
                        upper_bounds(i,j-1) = D_range(1);
                    end
                end
               
            %Transition at upper domain bound   
            elseif matching_position(2)
               lower_bounds(i,j) = transition_points(end);
               upper_bounds(i,j) = D_range(2); 
               if j < n_points
                        lower_bounds(i,j+1) = D_range(2);
                        upper_bounds(i,j+1) = D_range(2);
                end
            end
        end
        
    end
end

end

