function [dy,growth,pseudo_growth,total_Jg,osm] = pseudo_SS_RHS_ode15s(x)

%This script is the RHS of the pseudo-steady-state chemostat model used to
%compute the invasion fitness of a candidate strain.

global sspar

simu.rho = x(1:sspar.m);
simu.c = reshape(x(sspar.m+1:end),sspar.m+1,sspar.n_rxn+1);

%Define fundamental growth function
if sspar.no_saturation && ~sspar.no_thermo
    mm_fun = @(S,I,K,Keq) max([(sspar.nu).*(S - I./Keq),zeros(size(S))],[],2);
elseif sspar.no_saturation && sspar.no_thermo
    mm_fun = @(S,I,K,Keq) max([(sspar.nu).*(S),zeros(size(S))],[],2);
elseif ~sspar.no_saturation && sspar.no_thermo
    mm_fun = @(S,I,K,Keq) max([(sspar.nu/K).*(S)./(1 + S/K + I/K),zeros(size(S))],[],2);
else
    mm_fun = @(S,I,K,Keq) max([(sspar.nu/K).*(S - I./Keq)./(1 + S/K + I/K),zeros(size(S))],[],2);
end
    
%Extract intra and extracellular concentration
c_int = simu.c(1:end-1,:);
c_ext = simu.c(end,:);

%Compute matrix of growth rxn rates
Jg = zeros(sspar.m,sspar.n_rxn);
for i = 1:sspar.n_rxn
    Si = c_int(:,i);
    Ii = c_int(:,i+1);
    Jg(:,i) = mm_fun(Si,Ii,sspar.K(i),sspar.Keq(i));
end

%Multiply rate by enzyme levels
Jg = sspar.E.*Jg;
eff_Jg = Jg;

%Compute matrix of transport rates
Jt = sspar.beta.*sspar.T.*(repmat(c_ext,sspar.m,1) - c_int);

%Compute osmotic burdens
osm = sum(c_int,2);

%Compute cell growth rates, but make drho zero
total_Jg = sum(sspar.alpha.*eff_Jg,2);
total_Jg = max([total_Jg - sspar.min_growth,zeros(size(total_Jg))],[],2);
pseudo_growth = (1-sspar.osm_coeff*osm).*total_Jg;
growth = max([pseudo_growth,zeros(size(simu.rho))],[],2);
drho = zeros(size(simu.rho));

%Compute changes in concentrations
dc = zeros(size(simu.c));
dc(1:end-1,:) = Jt;

if ~sspar.neglect_pseudo_growth_dilution
dc(1:end-1,:) = dc(1:end-1,:) - repmat(growth,1,sspar.n_rxn+1).*simu.c(1:end-1,:);
end

dc(1:end-1,2:end) = dc(1:end-1,2:end) + Jg;

%Consumption in the intracellular compartment
dc(1:end-1,1:end-1) = dc(1:end-1,1:end-1) - Jg;

%Make final derivative vector
dy = [drho; dc(:)];

end

