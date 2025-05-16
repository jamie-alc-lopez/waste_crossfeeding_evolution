clear;clc

%This script takes in expressions for the optimal c0e of the P1C and P1P2
%consortia (computed by their respective "_lagrangian_derivative" scripts
%and makes a function of the ratio of the two for use in an SI fig

syms S alpha1 alpha2 beta D Vc Keq1 Keq2 gm positive

load('c0e_P1C.mat')
load('c0e_P1P2.mat')

P1C_P1P2_c0e_ratio = simplify(true_c0e_P1C/true_c0e_P1P2,'Steps',1000);

matlabFunction(P1C_P1P2_c0e_ratio,'File','P1C_P1P2_c0e_ratio',...
        'Vars',{S,alpha1,alpha2,beta,D,Keq1,Keq2,gm});