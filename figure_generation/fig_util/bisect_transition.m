function [current_bounds,transition_map] = bisect_transition(n_bisections,bfun,trn_id,transition_map,current_bounds)

%Loop through and perform bisection method
for j = 1:n_bisections
    
    %Get new middle value
    mid_D = mean(current_bounds(trn_id,:));
    mid_ESS_ind = bfun(mid_D);
    ind_triple = [transition_map(trn_id,1),mid_ESS_ind,transition_map(trn_id,2)];
    current_diff = abs(diff(ind_triple));
    
    %If there are are two transitions, add a secondary transition for
    %round two
    if sum(current_diff > 0) == 2
        disp('Multitransition detected, adding additional round of bisections.')
        
        %Note later transition for later round
        transition_map(end+1,1) = mid_ESS_ind;
        transition_map(end,2) = transition_map(trn_id,2);
        current_bounds(end+1,1) = mid_D;
        current_bounds(end,2) = current_bounds(trn_id,2);
        
        %Continue bisecting earlier transition
        current_bounds(trn_id,2) = mid_D;
        transition_map(trn_id,2) = mid_ESS_ind;
        
    %If only one transition found, pick either first or second interval
    elseif sum(current_diff > 0) == 1
        if current_diff(1) > 0
            current_bounds(trn_id,2) = mid_D;
        elseif current_diff(2) > 0
            current_bounds(trn_id,1) = mid_D;
        end
        
    %Flag warning if transition somehow does not appear in either half
    else
        disp('No transition detected within bisection.')
    end
    
end

end

