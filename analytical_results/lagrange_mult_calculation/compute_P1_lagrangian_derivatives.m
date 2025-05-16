clear;clc

%This script outputs the derivative of the Lagrangian for the P1
%optimization problem. From this point, we solved the system manually using
%the approach outlined in the SI

syms S alpha1 beta D rV Keq1 
syms E1P1 T0P1 T1P1 c0P1 c1P1 c0e c1e rhoP1
syms lambda

load('../simplified_model_symbolic_solution/P1_sol.mat')

h = -sym(1) + E1P1 + T0P1 + T1P1;

c0e_ss = sol.c0e(nsi);

dL_dE1P1 = simplify(diff(c0e_ss,E1P1)) - diff(lambda*h,E1P1)
dL_dT0P1 = simplify(diff(c0e_ss,T0P1)) - diff(lambda*h,T0P1)
dL_dT1P1 = simplify(diff(c0e_ss,T1P1)) - diff(lambda*h,T1P1)
dL_dlambda = simplify(diff(c0e_ss,lambda)) - diff(lambda*h,lambda)


