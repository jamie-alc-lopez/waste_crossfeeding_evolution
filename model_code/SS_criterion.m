function [value, isterminal, direction] = SS_criterion(t,y)
global par

%The function is the steady-state criterion for the chemostat simulation

dy = crossfeeding_RHS_ode15s(y);

%Check whether solution has reached tolerance, only looking at populations
%above the extinction threshold
drho = dy(1:par.m);
rho = y(1:par.m);
non_extinct_ind = rho >= par.extol; 

%Check if all species have gone extinct
if max(rho) < par.extol
    value = 0;
elseif datetime('now') > (par.start_time + minutes(par.max_runtime))
    value = 0;
    par.exceeded_time = 1;
else
    value = nanmax(abs(drho(non_extinct_ind)./rho(non_extinct_ind))) - par.tol;
end

isterminal = 1;   % Stop the integration
direction  = 0;

end