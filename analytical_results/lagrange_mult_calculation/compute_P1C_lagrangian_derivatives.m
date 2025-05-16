clear;clc

%This script outputs the derivative of the Lagrangian for the P1C
%optimization problem. From this point, we solved the system manually using
%the approach outlined in the SI

%This script also produces a function for the concentrations produced by
%the optimal P1C consortium

syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm positive
syms E1P1 T0P1 T1P1 E2C T1C T2C positive
syms c0P1 c1P1 c1C c2C c0e c1e c2e rhoP1 rhoC positive
syms lambdaP1 lambdaC

load('../simplified_model_symbolic_solution/P1C_sol.mat')

hP1 = -sym(1) + E1P1 + T0P1 + T1P1;
hC = -sym(1) + E2C + T1C + T2C;

c0e_ss = sol.c0e(nsi);

dL_dE1P1 = simplify(diff(c0e_ss,E1P1)) - diff(lambdaC*hC + lambdaP1*hP1,E1P1)
dL_dT0P1 = simplify(diff(c0e_ss,T0P1)) - diff(lambdaC*hC + lambdaP1*hP1,T0P1)
dL_dT1P1 = simplify(diff(c0e_ss,T1P1)) - diff(lambdaC*hC + lambdaP1*hP1,T1P1)
dL_dlambdaP1 = simplify(diff(c0e_ss,lambdaP1)) - diff(lambdaC*hC + lambdaP1*hP1,lambdaP1)

dL_dE2C = simplify(diff(c0e_ss,E2C)) - diff(lambdaC*hC + lambdaP1*hP1,E2C)
dL_dT1C = simplify(diff(c0e_ss,T1C)) - diff(lambdaC*hC + lambdaP1*hP1,T1C)
dL_dT2C = simplify(diff(c0e_ss,T2C)) - diff(lambdaC*hC + lambdaP1*hP1,T2C)
dL_dlambdaC = simplify(diff(c0e_ss,lambdaC)) - diff(lambdaC*hC + lambdaP1*hP1,lambdaC)

%% Plug in enzyme strategy solution to get steady-state concentrations

%Input manually computed optimal enzyme vals
E1P1_sol = 1/(1 + 1/sqrt(beta) + 1/sqrt(beta*Keq1));
T0P1_sol = 1/(1 + sqrt(beta) + 1/sqrt(Keq1));
T1P1_sol = 1/(1 + sqrt(Keq1) + sqrt(beta*Keq1));

E2C_sol = 1/(1 + 1/sqrt(beta) + 1/sqrt(beta*Keq2));
T1C_sol = 1/(1 + sqrt(beta) + 1/sqrt(Keq2));
T2C_sol = 1/(1 + sqrt(Keq2) + sqrt(beta*Keq2));

enz_vars = [E1P1,T0P1,T1P1,E2C,T1C,T2C];
enz_sols = [E1P1_sol,T0P1_sol,T1P1_sol,E2C_sol,T1C_sol,T2C_sol];

%Save just the c0e solution to make the c0e ratio function in a downstream
%script
true_c0e_P1C = simplify(subs(c0e_ss,enz_vars,enz_sols),'Steps',100);
save('c0e_P1C.mat','true_c0e_P1C')

% Function for internal values
conc_vec = [sol.c0e(nsi),sol.c1e(nsi),sol.c2e(nsi),sol.c0P1(nsi),sol.c1P1(nsi),sol.c1C(nsi),sol.c2C(nsi)];
true_conc_vec = simplify(subs(conc_vec,enz_vars,enz_sols));

matlabFunction(true_conc_vec,'File','optimal_P1C_concentrations',...
        'Vars',{S,alpha1,alpha2,beta,D,Keq1,Keq2,gm});

