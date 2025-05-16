clear;clc

%This script computes the analytical steady-state solution of the P1 only
%community in the simplified thermodynamic model with no internal dilution
%term. It produces a datafile containing the analytical expression and a
%MATLAB function that computes the steady state given parameter values.

syms S alpha1 beta D rV Keq1 gm positive real
syms E1P1 T0P1 T1P1 c0P1 c1P1 c0e c1e rhoP1 positive real

Jg1P1 = E1P1*(c0P1 - c1P1/Keq1);
gP1 = alpha1*(Jg1P1) - gm;
Jt0P1 = beta*T0P1*(c0P1 - c0e);
Jt1P1 = beta*T1P1*(c1P1 - c1e);
drho = (gP1 - D)*rhoP1;
dc0 = -Jg1P1 - Jt0P1;
dc1 = Jg1P1 - Jt1P1;
dc0e = S - D*c0e + rhoP1*Jt0P1*rV;
dc1e = -D*c1e + rhoP1*Jt1P1*rV;

vars = [c0P1 c1P1 c0e c1e rhoP1];
eqns = [drho == 0, dc0 == 0, dc1 == 0, dc0e == 0, dc1e == 0];

%symbolic solve
sol = solve(eqns,vars);

%Check for the P1 existence solution
nonzero_vars = {'rhoP1'};
nonzero_sol_ind = find_nonzero_sols(sol,nonzero_vars);
nsi = find(nonzero_sol_ind);

save('P1_sol.mat','sol','nsi');

%Generate function
eval_css = [simplify(sol.c0e(nsi)),simplify(sol.c1e(nsi)),...
    simplify(sol.rhoP1(nsi))];

matlabFunction(eval_css,'File','P1_equilibrium',...
        'Vars',{E1P1,T0P1,T1P1,S,alpha1,beta,D,rV,Keq1,gm});