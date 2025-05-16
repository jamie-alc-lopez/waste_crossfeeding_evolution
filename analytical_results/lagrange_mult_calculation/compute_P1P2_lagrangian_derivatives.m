clear;clc

%This script outputs the derivative of the Lagrangian for the P1P2
%optimization problem. From this point, we solved the system manually using
%the approach outlined in the SI

%This script also produces a function for the concentrations produced by
%the optimal P1P2 consortium

syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm positive
syms E1P1 T0P1 T1P1 E1P2 E2P2 T0P2 T2P2 positive 
syms c0P1 c1P1 c0P2 c1P2 c2P2 c0e c1e c2e rhoP1 rhoP2 positive
syms lambdaP1 lambdaP2

load('../simplified_model_symbolic_solution/P1P2_sol.mat')

hP1 = -sym(1) + E1P1 + T0P1 + T1P1;
hP2 = -sym(1) + E1P2 + E2P2 + T0P2 + T2P2;

c0e_ss = simplify(sol.c0e(nsi));

dL_dE1P1 = simplify(diff(c0e_ss,E1P1)) - diff(lambdaP1*hP1 + lambdaP2*hP2,E1P1)
dL_dT0P1 = simplify(diff(c0e_ss,T0P1)) - diff(lambdaP1*hP1 + lambdaP2*hP2,T0P1)
dL_dT1P1 = simplify(diff(c0e_ss,T1P1)) - diff(lambdaP1*hP1 + lambdaP2*hP2,T1P1)
dL_dlambdaP1 = simplify(diff(c0e_ss,lambdaP1)) - diff(lambdaP1*hP1 + lambdaP2*hP2,lambdaP1)

dL_dE1P2 = simplify(diff(c0e_ss,E1P2)) - diff(lambdaP1*hP1 + lambdaP2*hP2,E1P2)
dL_dE2P2 = simplify(diff(c0e_ss,E2P2)) - diff(lambdaP1*hP1 + lambdaP2*hP2,E2P2)
dL_dT0P2 = simplify(diff(c0e_ss,T0P2)) - diff(lambdaP1*hP1 + lambdaP2*hP2,T0P2)
dL_dT2P2 = simplify(diff(c0e_ss,T2P2)) - diff(lambdaP1*hP1 + lambdaP2*hP2,T2P2)
dL_dlambdaP2 = simplify(diff(c0e_ss,lambdaP2)) - diff(lambdaP1*hP1 + lambdaP2*hP2,lambdaP2)

%% Plug in enzyme strategy solution to get steady-state concentrations

%Manually input the optimal enzyme strategy for P1P2
E1P1_sol = 1/(1 + 1/sqrt(beta) + 1/sqrt(beta*Keq1));
T0P1_sol = 1/(1 + sqrt(beta) + 1/sqrt(Keq1));
T1P1_sol = 1/(1 + sqrt(Keq1) + sqrt(beta*Keq1));

E1P2_sol = 1/(1 + 1/sqrt(beta) + 1/sqrt(Keq1) + 1/(sqrt(Keq1*Keq2*beta)));
E2P2_sol = 1/(1 + sqrt(Keq1) + 1/sqrt(Keq2*beta) + sqrt(Keq1/beta));
T0P2_sol = 1/(1 + sqrt(beta) + 1/sqrt(Keq1*Keq2) + sqrt(beta/Keq1));
T2P2_sol = 1/(1 + sqrt(Keq1*Keq2) + sqrt(Keq2*beta) + sqrt(Keq1*Keq2*beta));

enz_vars = [E1P1,T0P1,T1P1,E1P2,E2P2,T0P2,T2P2];
enz_sols = [E1P1_sol,T0P1_sol,T1P1_sol,E1P2_sol,E2P2_sol,T0P2_sol,T2P2_sol];

%Save just the c0e solution to make the c0e ratio function in a downstream
%script
true_c0e_P1P2 = simplify(subs(c0e_ss,enz_vars,enz_sols),'Steps',100);
save('c0e_P1P2.mat','true_c0e_P1P2')

%Make a function for all concentration values
conc_vec = [sol.c0e(nsi),sol.c1e(nsi),sol.c2e(nsi),sol.c0P1(nsi),sol.c1P1(nsi),sol.c0P2(nsi),sol.c1P2(nsi),sol.c2P2(nsi)];
true_conc_vec = simplify(subs(conc_vec,enz_vars,enz_sols));
matlabFunction(true_conc_vec,'File','optimal_P1P2_concentrations',...
        'Vars',{S,alpha1,alpha2,beta,D,Keq1,Keq2,gm});

%%% Get analytical solution for population ratios
% 
% rhoP1 = sol.rhoP1(nsi);
% rhoP2 = sol.rhoP2(nsi);
% true_rhoP1_P1P2 = simplify(subs(rhoP1,enz_vars,enz_sols),'Steps',100);
% true_rhoP2_P1P2 = simplify(subs(rhoP2,enz_vars,enz_sols),'Steps',100);
% 
% c1 = sol.c1e(nsi);
% c2 = sol.c2e(nsi);
% 
% true_c1 = simplify(subs(c1,enz_vars,enz_sols),'Steps',100);
% true_c2 = simplify(subs(c2,enz_vars,enz_sols),'Steps',100);



