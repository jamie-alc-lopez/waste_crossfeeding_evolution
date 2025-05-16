clear;clc

%This script computes the analytical steady-state solution of the P2 only
%community in the simplified thermodynamic model with no internal dilution
%term. It produces a datafile containing the analytical expression and a
%MATLAB function that computes the steady state given parameter values.

%Set up system of equations
syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm
syms E1P2 E2P2 T0P2 T2P2 c0P2 c1P2 c2P2 c0e c2e rhoP2

Jg1P2 = E1P2*(c0P2 - c1P2/Keq1);
Jg2P2 = E2P2*(c1P2 - c2P2/Keq2);
gP2 = alpha1*Jg1P2 + alpha2*Jg2P2 - gm;
Jt0P2 = beta*T0P2*(c0P2 - c0e);
Jt2P2 = beta*T2P2*(c2P2 - c2e);
drhoP2 = (gP2 - D)*rhoP2;
dc0P2 = -Jg1P2 - Jt0P2;
dc1P2 = Jg1P2 - Jg2P2;
dc2P2 = Jg2P2 - Jt2P2;
dc0e = S - D*c0e + rhoP2*Jt0P2*rV;
dc2e = -D*c2e + rhoP2*Jt2P2*rV;

vars = [c0P2 c1P2 c2P2 c0e c2e rhoP2];
eqns = [drhoP2 == 0, dc0P2 == 0, dc1P2 == 0, dc2P2 == 0, dc0e == 0, dc2e == 0];

sol = solve(eqns,vars);

%Check for the P2 existence solution
nonzero_vars = {'rhoP2'};
nonzero_sol_ind = find_nonzero_sols(sol,nonzero_vars);
nsi = find(nonzero_sol_ind);

save('P2_sol.mat','sol','nsi');

%Generate function
eval_css = [simplify(sol.c0e(nsi)),simplify(sol.c2e(nsi)),...
    simplify(sol.rhoP2(nsi))];

matlabFunction(eval_css,'File','P2_equilibrium',...
        'Vars',{E1P2,E2P2,T0P2,T2P2,S,alpha1,alpha2,beta,D,rV,Keq1,Keq2,gm});