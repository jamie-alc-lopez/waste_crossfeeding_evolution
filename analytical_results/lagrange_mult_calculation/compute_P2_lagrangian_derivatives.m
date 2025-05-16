clear;clc

%This script outputs the derivative of the Lagrangian for the P1
%optimization problem. From this point, we solved the system manually using
%the approach outlined in the SI

syms S alpha1 beta D rV Keq1 Keq2 gm
syms E1P2 E2P2 T0P2 T2P2 c0P2 c1P2 c2P2 c0e c2e rhoP2
syms lambda

load('../simplified_model_symbolic_solution/P2_sol.mat')

h = -sym(1) + E1P2 + E2P2 + T0P2 + T2P2;

c0e_ss = sol.c0e(nsi);

dL_dE1P2 = simplify(diff(c0e_ss,E1P2)) - diff(lambda*h,E1P2)
dL_dE2P2 = simplify(diff(c0e_ss,E2P2)) - diff(lambda*h,E2P2)
dL_dT0P2 = simplify(diff(c0e_ss,T0P2)) - diff(lambda*h,T0P2)
dL_dT2P2 = simplify(diff(c0e_ss,T2P2)) - diff(lambda*h,T2P2)
dL_dlambda = simplify(diff(c0e_ss,lambda)) - diff(lambda*h,lambda)
