clear;clc

%This script computes the analytical steady-state solution of the P1P2 only
%community in the simplified thermodynamic model with no internal dilution
%term. It produces a datafile containing the analytical expression and a
%MATLAB function that computes the steady state given parameter values.

%Set up system of equations
syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm
syms E1P1 T0P1 T1P1 E1P2 E2P2 T0P2 T2P2 
syms c0P1 c1P1 c0P2 c1P2 c2P2 c0e c1e c2e rhoP1 rhoP2

Jg1P1 = E1P1*(c0P1 - c1P1/Keq1);
Jg1P2 = E1P2*(c0P2 - c1P2/Keq1);
Jg2P2 = E2P2*(c1P2 - c2P2/Keq2);
gP1 = alpha1*(Jg1P1) - gm;
gP2 = alpha1*Jg1P2 + alpha2*Jg2P2 - gm;

Jt0P1 = beta*T0P1*(c0P1 - c0e);
Jt1P1 = beta*T1P1*(c1P1 - c1e);
Jt0P2 = beta*T0P2*(c0P2 - c0e);
Jt2P2 = beta*T2P2*(c2P2 - c2e);

drhoP1 = (gP1 - D)*rhoP1;
drhoP2 = (gP2 - D)*rhoP2;

dc0P1 = -Jg1P1 - Jt0P1;
dc1P1 = Jg1P1 - Jt1P1;
dc0P2 = -Jg1P2 - Jt0P2;
dc1P2 = Jg1P2 -Jg2P2;
dc2P2 = Jg2P2 - Jt2P2;

dc0e = S - D*c0e + rhoP1*Jt0P1*rV + rhoP2*Jt0P2*rV;
dc1e = -D*c1e + rhoP1*Jt1P1*rV;
dc2e = -D*c2e + rhoP2*Jt2P2*rV;

vars = [rhoP1, rhoP2,c0P1,c1P1,c0P2,c1P2,c2P2,c0e,c1e,c2e];
eqns = [drhoP1 == 0, drhoP2== 0, dc0P1 == 0, dc1P1 == 0, dc0P2 == 0,...
    dc1P2 == 0, dc2P2 == 0, dc0e == 0, dc1e == 0, dc2e ==0];

%symbolic solve
sol = solve(eqns,vars);

%Check for the P1 and P2 coexistence solution
nonzero_vars = {'rhoP1','rhoP2'};
nonzero_sol_ind = find_nonzero_sols(sol,nonzero_vars);
nsi = find(nonzero_sol_ind);

save('P1P2_sol.mat','sol','nsi');

%Generate function
eval_css = [simplify(sol.c0e(nsi)),simplify(sol.c1e(nsi)),simplify(sol.c2e(nsi)),...
    simplify(sol.rhoP1(nsi)),simplify(sol.rhoP2(nsi))];

matlabFunction(eval_css,'File','P1P2_equilibrium',...
        'Vars',{E1P1,T0P1,T1P1,E1P2,E2P2,T0P2,T2P2,S,alpha1,alpha2,...
        beta,D,rV,Keq1,Keq2,gm});