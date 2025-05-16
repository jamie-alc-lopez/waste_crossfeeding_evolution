%% This script generates symbolic expression functions for the interior concentrations
% of a newly invading cell in the absence of growth. These functions are
% used in compute_pseudo_ss_analytical.m

%These functions differ from those used in the Langrange multiplier
%calculation as they are pseudo-equilibria, not the full equilibria

%% Compute versalitist sol
clear;clc
clear global
syms nu beta Keq1 Keq2
syms E1 E2 T0 T1 T2 c0 c1 c2 c0e c1e c2e

Jg1 = E1*(nu)*(c0 - c1/Keq1);
Jg2 = E2*(nu)*(c1 - c2/Keq2);
Jt0 = beta*T0*(c0 - c0e);
Jt1 = beta*T1*(c1 - c1e);
Jt2 = beta*T2*(c2 - c2e);
dc0 = -Jg1 - Jt0;
dc1 = Jg1 - Jg2 - Jt1;
dc2 = Jg2 - Jt2;

vars = [c0 c1 c2];
eqns = [dc0 == 0, dc1 == 0, dc2 == 0];

sol = solve(eqns,vars);

eval_css = [simplify(sol.c0),simplify(sol.c1),simplify(sol.c2)];

matlabFunction(eval_css,'File','vers_analytical_pseudo_equi',...
        'Vars',{E1,E2,T0,T1,T2,c0e,c1e,c2e,nu,beta,Keq1,Keq2});

%% Compute P1 sol
clear;clc
clear global
syms nu beta Keq1 Keq2
syms E1 E2 T0 T1 T2 c0 c1 c2 c0e c1e c2e

Jg1 = E1*(nu)*(c0 - c1/Keq1);
Jt0 = beta*T0*(c0 - c0e);
Jt1 = beta*T1*(c1 - c1e);
dc0 = -Jg1 - Jt0;
dc1 = Jg1 - Jt1;

vars = [c0 c1 c2];
eqns = [dc0 == 0, dc1 == 0];

sol = solve(eqns,vars);

eval_css = [simplify(sol.c0),simplify(sol.c1),simplify(sol.c2)];

matlabFunction(eval_css,'File','P1_analytical_pseudo_equi',...
    'Vars',{E1,E2,T0,T1,T2,c0e,c1e,c2e,nu,beta,Keq1,Keq2});

%% Compute P2 sol
clear;clc
clear global
syms nu beta Keq1 Keq2
syms E1 E2 T0 T1 T2 c0 c1 c2 c0e c1e c2e

Jg1 = E1*(nu)*(c0 - c1/Keq1);
Jg2 = E2*(nu)*(c1 - c2/Keq2);
Jt0 = beta*T0*(c0 - c0e);
Jt2 = beta*T2*(c2 - c2e);
dc0 = -Jg1 - Jt0;
dc1 = Jg1 - Jg2;
dc2 = Jg2 - Jt2;

vars = [c0 c1 c2];
eqns = [dc0 == 0, dc1 == 0, dc2 == 0];

sol = solve(eqns,vars);

eval_css = [simplify(sol.c0),simplify(sol.c1),simplify(sol.c2)];

matlabFunction(eval_css,'File','P2_analytical_pseudo_equi',...
        'Vars',{E1,E2,T0,T1,T2,c0e,c1e,c2e,nu,beta,Keq1,Keq2});

%% Compute C sol
clear;clc
clear global
syms nu beta Keq1 Keq2
syms E1 E2 T0 T1 T2 c0 c1 c2 c0e c1e c2e

Jg2 = E2*(nu)*(c1 - c2/Keq2);
Jt1 = beta*T1*(c1 - c1e);
Jt2 = beta*T2*(c2 - c2e);
dc1 = -Jg2 - Jt1;
dc2 = Jg2 - Jt2;

vars = [c0 c1 c2];
eqns = [dc1 == 0, dc2 == 0];

sol = solve(eqns,vars);

eval_css = [simplify(sol.c0),simplify(sol.c1),simplify(sol.c2)];

matlabFunction(eval_css,'File','C_analytical_pseudo_equi',...
        'Vars',{E1,E2,T0,T1,T2,c0e,c1e,c2e,nu,beta,Keq1,Keq2});

    

%% Test against numerical pseudo equilibrium
clear;clc
clear global
global par
load('../parameters/default_parameter_struct.mat');
par.no_saturation = 1;
par.no_thermo = 0;
par.no_T1 = 0;
%simu.c = [10 0.2 0.1];
%ET = [0.2, 0, 0.3, 0.2, 0.3];
simu.c = [16.1914780274016,3.81153171529499,2.25122047482084];
ET = [0.0852904436099165,0,0.0670933236382738,0.847616232751810,0];

par.osm_coeff = 0.1;
tic
[pseudo_growth1,growth1,c_ss1] = compute_pseudo_SS_ode15s(ET,simu);
toc

tic
[pseudo_growth2,growth2,c_ss2] = compute_pseudo_SS_analytical(ET,simu);
toc

tic
[pseudo_growth3,growth3,c_ss3] = compute_pseudo_SS_fsolve(ET,simu);
toc

