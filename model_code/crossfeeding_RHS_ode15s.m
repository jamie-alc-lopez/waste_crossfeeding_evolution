function dy = crossfeeding_RHS_ode15s(x)

%This script is the RHS of the full chemostat model

global par

%Reassemble the state variable struct
simu.rho = x(1:par.m);
simu.c = reshape(x(par.m+1:end),par.m+1,par.n_rxn+1);

%Define fundamental growth function
if par.no_saturation && ~par.no_thermo
    mm_fun = @(S,I,K,Keq) max([(par.nu).*(S - I./Keq),zeros(size(S))],[],2);
elseif par.no_saturation && par.no_thermo
    mm_fun = @(S,I,K,Keq) max([(par.nu).*(S),zeros(size(S))],[],2);
elseif ~par.no_saturation && par.no_thermo
    mm_fun = @(S,I,K,Keq) max([(par.nu/K).*(S)./(1 + S/K + I/K),zeros(size(S))],[],2);
else
    mm_fun = @(S,I,K,Keq) max([(par.nu/K).*(S - I./Keq)./(1 + S/K + I/K),zeros(size(S))],[],2);
end

%Extract intra and extracellular concentration
c_int = simu.c(1:end-1,:);
c_ext = simu.c(end,:);

%Compute matrix of growth rxn rates
Jg = zeros(par.m,par.n_rxn);
for i = 1:par.n_rxn
    Si = c_int(:,i);
    Ii = c_int(:,i+1);
    Jg(:,i) = mm_fun(Si,Ii,par.K(i),par.Keq(i));
end


%Multiply rate by enzyme levels
Jg = par.E.*Jg;
eff_Jg = Jg;

%Compute matrix of transport rates
Jt = par.beta.*par.T.*(repmat(c_ext,par.m,1) - c_int);

%Compute osmotic burdens
osm = sum(c_int,2);

%Compute cell growth rates, using minimum growth and osmotic penalty
total_Jg = sum(par.alpha.*eff_Jg,2);
total_Jg = max([total_Jg - par.min_growth,zeros(size(total_Jg))],[],2);
pseudo_growth =(1-par.osm_coeff*osm).*total_Jg;
growth = max([pseudo_growth,zeros(size(simu.rho))],[],2);
drho = growth.*simu.rho;

%Compute changes in concentrations (note that the osmolarity impacts yield only)
dc = zeros(size(simu.c));

%Transport - dilution in intracellular compartment
dc(1:end-1,:) = Jt - repmat(growth,1,par.n_rxn+1).*simu.c(1:end-1,:);

%Production in intracellular compartment
dc(1:end-1,2:end) = dc(1:end-1,2:end) + Jg;

%Transport in extracellular compartment
dc(end,:) = -sum(par.VcVC*Jt.*repmat(simu.rho,1,size(simu.c,2)),1);

%Consumption in the intracellular compartment
dc(1:end-1,1:end-1) = dc(1:end-1,1:end-1) - Jg;

%Add in inputs/outputs
drho = drho - simu.rho*par.D;
dc(end,:) = dc(end,:) - c_ext*par.D;
dc(end,1) = dc(end,1) + par.supply;

if isfield(par,'supply2')
    dc(end,2) = dc(end,2) + par.supply2;
end

if isfield(par,'supply3')
    dc(end,3) = dc(end,3) + par.supply3;
end
    
%Make final derivative vector
dy = [drho; dc(:)];

end

