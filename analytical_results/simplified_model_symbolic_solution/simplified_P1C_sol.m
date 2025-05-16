clear;clc

%This script computes the analytical steady-state solution of the P1+C
%community in the simplified thermodynamic model with no internal dilution
%term. It produces a datafile containing the analytical expression and a
%MATLAB function that computes the steady state given parameter values.

%Set system of equations
syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm
syms E1P1 T0P1 T1P1 E2C T1C T2C 
syms c0P1 c1P1 c1C c2C c0e c1e c2e rhoP1 rhoC 

Jg1P1 = E1P1*(c0P1 - c1P1/Keq1);
Jg2C = E2C*(c1C - c2C/Keq2);
gP1 = alpha1*(Jg1P1) - gm;
gC = alpha2*(Jg2C) - gm;
Jt0P1 = beta*T0P1*(c0P1 - c0e);
Jt1P1 = beta*T1P1*(c1P1 - c1e);
Jt1C = beta*T1C*(c1C - c1e);
Jt2C = beta*T2C*(c2C - c2e);
drhoP1 = (gP1 - D)*rhoP1;
drhoC = (gC - D)*rhoC;
dc0P1 = -Jg1P1 - Jt0P1;
dc1P1 = Jg1P1 - Jt1P1;
dc1C = -Jg2C - Jt1C;
dc2C = Jg2C - Jt2C;
dc0e = S - D*c0e + rhoP1*Jt0P1*rV;
dc1e = -D*c1e + rhoP1*Jt1P1*rV + rhoC*Jt1C*rV;
dc2e = -D*c2e + rhoC*Jt2C*rV;

vars = [rhoP1, rhoC,c0P1,c1P1,c1C,c2C,c0e,c1e,c2e];
eqns = [drhoP1 == 0, drhoC== 0, dc0P1 == 0, dc1P1 == 0, dc1C == 0, ...
    dc2C == 0, dc0e == 0, dc1e == 0, dc2e ==0];

%symbolic solve
sol = solve(eqns,vars,'MaxDegree', 4);

%Check for the P1 and C coexistence solution
nonzero_vars = {'rhoP1','rhoC'};
nonzero_sol_ind = find_nonzero_sols(sol,nonzero_vars);
nsi = find(nonzero_sol_ind);

save('P1C_sol.mat','sol','nsi');

%Generate function
eval_css = [simplify(sol.c0e(nsi)),simplify(sol.c1e(nsi)),simplify(sol.c2e(nsi)),...
    simplify(sol.rhoP1(nsi)),simplify(sol.rhoC(nsi))];

matlabFunction(eval_css,'File','P1C_equilibrium',...
        'Vars',{E1P1,T0P1,T1P1,E2C,T1C,T2C,S,alpha1,alpha2,beta,D,rV,Keq1,Keq2,gm});