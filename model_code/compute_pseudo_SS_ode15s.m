function [pseudo_growth,growth,c_ss,total_Jg,osm] = compute_pseudo_SS_ode15s(ET,simu)

%This script computes the pseudo-steady-state of an invading strain. It is
%used during the invading strain optimization procedure.

global sspar
global par

%Define minimal simulation under a separate parameter structure
sspar = par;
ET(ET < 0) = 0;

if isfield(sspar,'only_P1')
    if sspar.only_P1
        ET = ET.*[1 0 1 1 0];
    end
end
if isfield(sspar,'only_P2')
    if sspar.only_P2
        ET = ET.*[1 1 1 0 1];
    end
end
if isfield(sspar,'only_C')
    if sspar.only_C
        ET = ET.*[0 1 0 1 1];
    end
end

%Initiate simulation
sspar.E = ET(1:sspar.n_rxn);
sspar.T = ET(sspar.n_rxn+1:end);
sspar.m = 1;
simu.rho = 1;

%Use the analytical prediction as an initial conditions
c_new = zeros(size(simu.c(end,:)));
%[~,~,c_new] = compute_pseudo_SS_analytical(ET,simu);
%c_new(isnan(c_new)) = 0;
simu.c = [c_new;simu.c(end,:)];

%Run a minimal simulation until pseudo-equilibrium is reached
y0 = [simu.rho;simu.c(:)];
Opt = odeset('NonNegative',1:length(y0));

tspan = [0,1e5];

odefun = @(t,y) pseudo_SS_RHS_ode15s(y);

loop_var = true;
tn = 0;
while loop_var
    %Compute RHS
    dy = pseudo_SS_RHS_ode15s(y0);
    int_conc_ind = false(size(y0));

    %Select only meaningful variables to compare to tolerance
    int_conc_ind(2:2:end) = true;
    tol_conc_ind = int_conc_ind & (y0 > 1e-3);

    %Compute max dy
    %max_dy = nanmax(abs(dy(tol_conc_ind)./y0(tol_conc_ind)));
    max_dy = nanmax(abs(dy));

    if max_dy < sspar.pseudo_tol
        yfinal = y0;
        loop_var = false;
    elseif tn > sspar.max_pseudo_rounds
        yfinal = y0;
        loop_var = false;
    else
        [~,y] = ode15s(odefun,tspan,y0,Opt);
        y0 = y(end,:)';
    end
    tn = tn + 1;
end
c_ss = reshape(yfinal(sspar.m+1:end),sspar.m+1,sspar.n_rxn+1);
c_ss = c_ss(1,:);
[~,growth,pseudo_growth,total_Jg,osm] = odefun(0,yfinal);

end

