function [optimal_ET,opt_growth,c_ss,optimal_ET_cell,opt_growth_vec] = compute_optimal_strain(simu)

%This script computes the optimal strain to invade a given environment.

global par

%Apply defaults if not specified
if ~isfield(par,'no_saturation')
    par.no_saturation = false;
end
if ~isfield(par,'no_thermo')
    par.no_thermo = false;
end
if ~isfield(par,'fast_pseudo')
    par.fast_pseudo = false;
end

%Make -fitness the objective
if par.fast_pseudo
    fun = @(x) -compute_pseudo_SS_analytical(x,simu);
else
    fun = @(x) -compute_pseudo_SS_ode15s(x,simu);
end

%Set-up optimizer parameters
A = [];
b = [];
Aeq = ones(1,2*par.n_rxn + 1);
beq = 1;
%lb = zeros(1,2*par.n_rxn + 1);

%Set upper bounds for enzymes


ub_mat = [1 1 1 0 1; 0 1 0 1 1; 1 0 1 1 0];
%ub_mat = [1 1 1 1 1; 1 1 1 0 1; 0 1 0 1 1; 1 0 1 1 0];

lb_mat = ub_mat*0; 

%Get strategy subset if only certain consortia allowed
if isfield(par,'only_P1P2')
    if par.only_P1P2
        strat_subset = [1;3];
    end
end

if isfield(par,'only_P1C')
    if par.only_P1C
        strat_subset = [2;3];
    end
end

%Strategy subset
if exist('strat_subset')
    ub_mat = ub_mat(strat_subset,:);
end


%List of initial conditions and algorithms
x0_mat = [1/5, 1/5, 1/5, 1/5, 1/5; ...
    1/3-2e-15, 1e-15, 1/3, 1/3, 1e-15;...
    1e-15, 1/3, 1e-15, 1/3, 1/3;...
    1/4, 1/4, 1/4, 1e-15, 1/4;...
    1/2-1/300, 1e-15, 1/300, 1/2-2e-15, 1e-15;...
    1/4-1e-15, 1/4, 1/300, 1e-15, 1/2-1/300;...
    1e-15, 1/3-2e-15, 1e-15, 1/300, 2/3-1/300;...
    0.105750129058066,0,0.0313763130798441,0.862873557862090,0;...
    0,0.0798750122534948,0,0.402613591877500,0.517511395869006;...
    0.8751, 0, 0.0624370206199267, 0.0624370115315309, 0;...
    0.94-0.001, 0.001, 0.04, 0.02, 0];

%Generate linear analytical predictions
E1P1_sol = 1/(1 + 1/sqrt(par.beta) + 1/sqrt(par.beta*par.Keq(1)));
T0P1_sol = 1/(1 + sqrt(par.beta) + 1/sqrt(par.Keq(1)));
T1P1_sol = 1/(1 + sqrt(par.Keq(1)) + sqrt(par.beta*par.Keq(1)));
E1P2_sol = 1/(1 + 1/sqrt(par.beta) + 1/sqrt(par.Keq(1)) +...
    1/(sqrt(par.Keq(1)*par.Keq(2)*par.beta)));
E2P2_sol = 1/(1 + sqrt(par.Keq(1)) + 1/sqrt(par.Keq(2)*par.beta) + ...
    sqrt(par.Keq(1)/par.beta));
T0P2_sol = 1/(1 + sqrt(par.beta) + 1/sqrt(par.Keq(1)*par.Keq(2)) + ...
    sqrt(par.beta/par.Keq(1)));
T2P2_sol = 1/(1 + sqrt(par.Keq(1)*par.Keq(2)) + ...
    sqrt(par.Keq(2)*par.beta) + sqrt(par.Keq(1)*par.Keq(2)*par.beta));

linear_P1 = [E1P1_sol 0 T0P1_sol T1P1_sol 0];
linear_C = [0 E1P1_sol 0 T0P1_sol T1P1_sol];
linear_P2 = [E1P2_sol E2P2_sol T0P2_sol 0 T2P2_sol];
x0_mat = [x0_mat; linear_P1; linear_C; linear_P2];

x0_mat = [x0_mat; [par.E par.T]];

%Set of algorithms
algorithms = {'interior-point','sqp'};

%Perform constrained optimization across different initial conditions
options = optimoptions(@fmincon);
options.Display =  par.reporting;
options.OptimalityTolerance = 1e-8;
iter = 1;
for i = 1:size(x0_mat,1)
    for j = 1:length(algorithms)
        options.Algorithm = algorithms{j};
        for k = 1:size(ub_mat,1)
            algo_cell{iter} = algorithms{j};
            IC_vec(iter) = i;
            ub_vec(iter) = k;
            [optimal_ET_cell{iter},opt_growth_vec(iter)] = ...
                fmincon(fun,x0_mat(i,:),A,b,Aeq,beq,lb_mat(k,:),ub_mat(k,:),[],options);
            iter = iter + 1;
        end
    end
end

%Find the best strain from this process
max_ind = find(opt_growth_vec == min(opt_growth_vec));
optimal_ET = optimal_ET_cell{max_ind(1)};
opt_growth = opt_growth_vec(max_ind(1));

%get ss values of c
if par.fast_pseudo
    [~,~,c_ss] = compute_pseudo_SS_analytical(optimal_ET,simu);
else
    [~,~,c_ss] = compute_pseudo_SS_ode15s(optimal_ET,simu);
end

%Plot optimization results if requested
if par.plt_opt_results
    figure
    violinplot(-opt_growth_vec,ub_vec)
    ylabel('Invasion growth rate')
    xticklabels({'Only P2','Only C','Only P1'})
    set(gca,'FontSize',15)
end

end

