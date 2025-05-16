clear;clc

%This script uses the evolutionarily stable strategies from the Lagrange
%multiplier calculation to compute the expression for the difference 
%between the optimal c0e_P1C and c0e_P1P2 

syms S alpha1 alpha2 beta D rV Keq1 Keq2 gm real positive
syms alpha Keq real positive

%Load c0e solutions calculated from Lagrange multiplier derivation
P1P2_sol_data = load('c0e_P1P2.mat');
P1C_sol_data = load('c0e_P1C.mat');

%Full difference expression
P1C_P1P2_c0e_diff = ... 
    simplify(P1C_sol_data.true_c0e_P1C - P1P2_sol_data.true_c0e_P1P2,...
    'steps',5000);

%Simplified expression with equal Keq and alpha
simplified_c0e_diff = simplify(subs(P1C_P1P2_c0e_diff,...
    [Keq1,Keq2,alpha1,alpha2],[Keq, Keq, alpha, alpha]),'steps',5000);

%Test whether simplified expression has region where difference is negative
equal_Keq_alpha_always_positive = isAlways(simplified_c0e_diff > 0);

%Simplified expression with D + gm = 1 for numerical evaluation
numerical_P1P2_c0e_diff = subs(P1C_P1P2_c0e_diff,...
    [D,gm],[sym(0.5),sym(0.5)]);
matlabFunction(numerical_P1P2_c0e_diff,'File','numerical_P1P2_c0e_diff',...
        'Vars',{alpha1,alpha2,Keq1,Keq2,beta});




%% Numerical evaluate whether there is a case where the difference is negative

clear;clc

%Define parameter range to test
n = 100;
Keq1 = logspace(-5,5,n);
Keq2 = logspace(-5,5,n);
alpha1 = logspace(-5,5,n);
alpha2 = logspace(-5,5,n);
beta = logspace(-5,5,n);

%Make a multidimensional grid of all parameters except Keq1
[alpha1_mat,alpha2_mat,Keq2_mat,beta_mat] = ...
    ndgrid(alpha1,alpha2,Keq2,beta);

%Subdivide the calculation to be more memory efficient
for i = 1:length(Keq1)

    %Select ith Keq1 value
    Keq1_i = Keq1(i);

    %Compute all c0e difference values with ith Keq1
    temp_diff_mat = numerical_P1P2_c0e_diff(alpha1_mat,alpha2_mat,Keq1_i,Keq2_mat,beta_mat);

    %Find minimum difference of this calculation round
    partial_min_diff(i) = min(temp_diff_mat(:));

    disp(['Completed calculation round ',num2str(i),' of ',num2str(n),'.'])
end

min_diff = min(partial_min_diff);