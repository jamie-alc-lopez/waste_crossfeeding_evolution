function [ESS_ind,ESS_type,optimal_ET,optimal_c1e,ss_vals,optimal_c1e_vec,xopt_cell] = ...
    compute_ESS_consortium(S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm)

%This function finds the evolutionarily stable consortium based on the
%thermodynamic toxicity model neglecting growth dilution. 

%Define possible consortia
possible_types = {'P1','P2','P1P2','P1C'};

%Define optimization function for possible consortia
fun_cell{1} = @(x) P1_equilibrium(x(1),x(2),x(3),S,alpha1,beta,D,Vc,Keq1,gm);
Aeq_cell{1} = ones(1,3);

fun_cell{2} = @(x) P2_equilibrium(x(1),x(2),x(3),x(4),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
Aeq_cell{2} = ones(1,4);

fun_cell{3} = @(x) P1P2_equilibrium(x(1),x(2),x(3),x(4),x(5),x(6),x(7),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
Aeq_cell{3} = [1 1 1 0 0 0 0; 0 0 0 1 1 1 1];

fun_cell{4} = @(x) P1C_equilibrium(x(1),x(2),x(3),x(4),x(5),x(6),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
Aeq_cell{4} = [1 1 1 0 0 0; 0 0 0 1 1 1];

%Find optimal of each consortia
xopt_cell = cell(size(fun_cell));
optimal_c1e_vec = zeros(size(fun_cell));
for i = 1:length(fun_cell)
    %Create function returning c1e as scalar
    fun = fun_cell{i};
    Aeq = Aeq_cell{i};
    minfun = @(x) subsref(fun(x),struct('type','()','subs',{{1}}));
    A = [];
    b = [];
    beq = ones(size(Aeq,1),1);
    lb = zeros(1,size(Aeq,2));
    ub = ones(1,size(Aeq,2));
    x0 = ub'*0.5;
    [xopt_cell{i},optimal_c1e_vec(i)] = fmincon(minfun,x0,A,b,Aeq,beq,lb,ub);
end

%Identify ESS
ESS_ind = find(optimal_c1e_vec == min(optimal_c1e_vec));
ESS_type = possible_types{ESS_ind};
xopt = xopt_cell{ESS_ind};
optimal_c1e = optimal_c1e_vec(ESS_ind);

%Format enzyme matrix and get biomass

if ESS_ind == 1
    optimal_ET = [xopt(1) 0 xopt(2) xopt(3) 0];
    ss_vals = P1_equilibrium(xopt(1),xopt(2),xopt(3),S,alpha1,beta,D,Vc,Keq1,gm);
elseif ESS_ind == 2
    optimal_ET = [xopt(1) xopt(2) xopt(3) 0 xopt(4)];
    ss_vals = P2_equilibrium(xopt(1),xopt(2),xopt(3),xopt(4),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
elseif ESS_ind == 3
    optimal_ET = [xopt(1) 0 xopt(2) xopt(3) 0; xopt(4) xopt(5) xopt(6) 0 xopt(7)];
    ss_vals = P1P2_equilibrium(xopt(1),xopt(2),xopt(3),xopt(4),xopt(5),xopt(6),xopt(7),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
elseif ESS_ind == 4
    optimal_ET = [xopt(1) 0 xopt(2) xopt(3) 0; 0 xopt(4) 0 xopt(5) xopt(6)];
    ss_vals = P1C_equilibrium(xopt(1),xopt(2),xopt(3),xopt(4),xopt(5),xopt(6),S,alpha1,alpha2,beta,D,Vc,Keq1,Keq2,gm);
end

if sum(ss_vals < 0) > 0
    ESS_ind = 0;
    ESS_type = 'ext';
end

end